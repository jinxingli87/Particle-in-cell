!-----------------------------------------------------------------------
! PIC simulation of the Generation of banded chorus waves 
! Wirtten by Jinxing Li, UCLA at 6/9/2017
! Based on PIC code composed by Viktor K. Decyk, UCLA
! Skeleton 1-2/2D Darwin PIC code
program dpic1
    implicit none
    real,external:: ranorm
    ! indx = exponent which determines grid points in x direction:
    ! nx = 2**indx.
    integer, parameter :: indx =   10
    ! npx = number of electrons distributed in x direction.
    integer, parameter :: npx =  2**24
    ! tend = time at end of simulation, in units of plasma frequency.
    ! dt = time interval between successive calculations.
    ! qme = charge on electron, in units of e.
    real, parameter :: tend = 40000.0, dt = 0.2, qme = -1.0
    !================================================================

    integer,parameter :: nkappa=2
    real, dimension(nkappa) :: vtx,vty,vx0,kappa
    integer :: ne(nkappa)
    integer,parameter :: T_Interval=4
    integer,parameter :: Np_Interval=16
    !================================================================
    ! ax = smoothed particle size in x direction
    ! ci = reciprocal of velocity of light.
    real :: ax = .912871, ci = 0.05
    ! idimp = number of particle coordinates = 4
    ! ipbc = particle boundary condition: 1 = periodic
    ! sortime = number of time steps between standard electron sorting
    integer :: idimp = 4, ipbc = 1, sortime = 50
    ! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
    real :: omx = 0.193, omy = 0.0518, omz = 0.0 !om=0.2; wna=15 deg
    ! ndc = number of corrections in darwin iteration
    integer :: ndc = 1
    ! wke/we = particle kinetic/electrostatic field energy
    ! wf/wm/wt = magnetic field/transverse electric field/total energy
    real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
    real :: zero = 0.0
    ! declare scalars for standard code
    integer :: k
    integer :: np, nx, nxh, nxe, nxeh
    integer :: nx1, ntime, nloop, isign
    real :: qbme, affp, q2m0, wpm, wpmax, wpmin
    !
    ! declare arrays for standard code:
    ! part, part2 = particle arrays
    real, dimension(:,:), pointer :: part, part2, tpart
    ! qe = electron charge density with guard cells
    ! fxe = smoothed longitudinal electric field with guard cells
    real, dimension(:), pointer :: qe, fxe
    ! cue = electron current density with guard cells
    ! dcu = acceleration density with guard cells
    ! cus = transverse electric field with guard cells
    ! amu = momentum flux with guard cells
    real, dimension(:,:), pointer :: cue, dcu, cus, amu
    ! exyze/byze = smoothed electric/magnetic field with guard cells
    real, dimension(:,:), pointer :: exyze, byze
    ! ffc, ffe = form factor arrays for poisson solvers
    complex, dimension(:), pointer :: ffc, ffe
    ! mixup = bit reverse table for FFT
    integer, dimension(:), pointer :: mixup
    ! sct = sine/cosine table for FFT
    complex, dimension(:), pointer :: sct
    ! npic = scratch array for reordering particles
    integer, dimension(:), pointer :: npic
    ! gxe, gyze = scratch arrays for fft
    real, dimension(:), pointer :: gxe
    real, dimension(:,:), pointer :: gyze
    !
    ! declare and initialize timing data
    real :: time
    integer, dimension(4) :: itime
    real :: tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0
    real :: tdjpost = 0.0, tdcjpost = 0.0, tpush = 0.0, tsort = 0.0
    double precision :: dtime

    ! =====================Added by Jinxing Li==================
    integer :: kk,i,j
    integer :: nloop_percent
    character filename *50,dir1 *50,filename1 *50, filename2 *50
    character filename3 *50,filename4 *50
    logical :: dir_e
    integer :: modesx, modesxd
    complex, dimension(:), pointer :: potc, pott
    complex, dimension(:,:), pointer :: vpotc, vpott,Bw_k
    complex, dimension(:), pointer :: Ex_k

    integer,parameter:: nv=91,nth=16
    !nth must be even number, so the number of pitch angle sectors,
    !nth-1, is odd, and there is a secotor covering 90 degree.
    real:: vArr(nv),thArr(nth),vArr1(nv-1),thArr1(nth-1)
    real:: npsd(nv-1,nth-1),fpsd(nv-1,nth-1),fvpara(nv-1),fvperp(nv-1)
    real log_vmin,log_vmax

    dir1='/home/jli/results/pic/whistler_kappa39/'
    inquire(file=trim(dir1),exist=dir_e)
    if(dir_e) then
    else
        call system('mkdir '//trim(dir1))
    end if

    write(filename1,*),trim(dir1),'pot.dat'
    write(filename2,*),trim(dir1),'para.txt'
    write(filename3,*),trim(dir1),'field.dat'
    write(filename4,*),trim(dir1),'edist.txt'
    open(97,file=trim(adjustl(filename1)),status='replace',&
        &action='write',form='unformatted',access='stream')
    open(100,file=trim(adjustl(filename3)),status='replace',&
        &action='write',form='unformatted',access='stream')
    open(150,file=trim(adjustl(filename4)))
    ! ==========================================================
    !
    ! initialize scalars for standard code
    ! np = total number of particles in simulation
    ! nx = number of grid points in x direction
    np = npx; nx = 2**indx; nxh = nx/2
    nxe = nx + 2; nxeh = nxe/2; nx1 = nx + 1
    ! nloop = number of time steps in simulation
    ! ntime = current time step
    nloop = tend/dt + .0001; ntime = 0
    qbme = qme
    affp = real(nx)/real(np)
    !
    ! allocate data for standard code
    allocate(part(idimp,np))
    if (sortime > 0) allocate(part2(idimp,np))
    allocate(qe(nxe),fxe(nxe))
    allocate(cue(2,nxe),dcu(2,nxe),cus(2,nxe),amu(2,nxe))
    allocate(exyze(3,nxe),byze(2,nxe))
    allocate(ffc(nxh),ffe(nxh),mixup(nxh),sct(nxh))
    allocate(npic(nx1))
    allocate(gxe(nxe),gyze(2,nxe))
    !
    ! prepare fft tables
    call WFFT1RINIT(mixup,sct,indx,nxh)
    ! calculate form factor: ffc
    isign = 0
    call POIS1(qe,fxe,isign,ffc,ax,affp,we,nx)
    !===========================================================
    ! initialize electrons
    ne(1)=int(np*0.8)
    ne(2)=np-ne(1)
    vtx=[0.20,0.4]
    vty=[0.20,0.8]!this is in either of the perpendicular direction
    vx0=[0.0,0.0]
    kappa=[4.0,1.5]

    log_vmin=-2.0
    log_vmax=1.3
    do i=1,nv
        vArr(i)=10**((log_vmax-log_vmin)/(nv-1)*real(i-1)+log_vmin)
    end do
    do i=1,nth
        thArr(i)=3.1415926/(nth-1)*(i-1)
    end do
    !call DISTR1H(part,vtx,vty,vtz,vx0,vy0,vz0,npx,idimp,np,nx,ipbc)
    call edist(part,nkappa,nx,np,ne,vtx,vty,vx0,kappa,log_vmin-1,log_vmax,nv*2-1,omx,omy)
    !In the initial setting, the velocity lower limit should be down 1 magnitude
    !and the nv doubles, so that the distribution near PA=0, 90 deg have
    !sufficient particles

    !===========================================================
    !
    ! find maximum and minimum initial electron density
    qe = 0.0
    call GPOST1L(part,qe,qme,np,idimp,nxe)
    call AGUARD1L(qe,nx,nxe)
    call FWPMINMX1(qe,qbme,wpmax,wpmin,nx,nxe)
    wpm = 0.5*(wpmax + wpmin)*affp
    ! accelerate convergence: update wpm
    if (wpm <= 10.0) wpm = 0.75*wpm
    write (*,*) 'wpm=',wpm
    q2m0 = wpm/affp
    ! calculate form factor: ffe
    isign = 0
    call EPOIS13(dcu,cus,isign,ffe,ax,affp,wpm,ci,wf,nx,nxeh,nxh)
    !
    ! initialize transverse electric field
    cus = 0.0
    !
    ! * * * start main iteration loop * * *
    !
    !==========    Edited by Jinxing Li   =======================
    modesx = nx/4; modesxd = modesx
    allocate(potc(nxeh),vpotc(2,nxeh))
    allocate(pott(modesxd),vpott(2,modesxd),Bw_k(2,nxeh))
    allocate(Ex_k(nxeh))
    !============================================================
    !500 if (nloop <= ntime) go to 2000
    nloop_percent=int(nloop/100)
    if (nloop_percent <1) then
        nloop_percent=1
    end if
    do ntime=0,nloop-1
        if (mod(ntime,nloop_percent) == 0) then
            print *, int(ntime/nloop_percent), " % completed"
        end if

        !     write (*,*) 'ntime = ', ntime
        !
        ! deposit current with standard procedure: updates cue
        call dtimer(dtime,itime,-1)
        cue = 0.0
        call GJPOST1L(part,cue,qme,zero,np,idimp,nx,nxe,ipbc)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tdjpost = tdjpost + time
        !
        ! deposit charge with standard procedure: updates qe
        call dtimer(dtime,itime,-1)
        qe = 0.0
        call GPOST1L(part,qe,qme,np,idimp,nxe)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tdpost = tdpost + time
        !
        ! add guard cells with standard procedure: updates qe, cue
        call dtimer(dtime,itime,-1)
        call AGUARD1L(qe,nx,nxe)
        call ACGUARD1L(cue,nx,nxe)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tguard = tguard + time
        !
        ! transform charge to fourier space with standard procedure:
        ! updates qe, gxe
        call dtimer(dtime,itime,-1)
        isign = -1
        call FFT1RXX(qe,gxe,isign,mixup,sct,indx,nxe,nxh)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tfft = tfft + time
        !
        ! calculate longitudinal force/charge in fourier space with standard
        ! procedure: updates fxe, we
        call dtimer(dtime,itime,-1)
        isign = -1
        call POIS1(qe,fxe,isign,ffc,ax,affp,we,nx)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tfield = tfield + time
        !
        ! transform longitudinal electric force to real space with standard
        ! procedure: updates fxe, gxe
        call dtimer(dtime,itime,-1)
        isign = 1
        call FFT1RXX(fxe,gxe,isign,mixup,sct,indx,nxe,nxh)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tfft = tfft + time
        !
        ! transform current to fourier space with standard procedure:
        ! updates cue, gyze
        call dtimer(dtime,itime,-1)
        isign = -1
        call FFT1R2X(cue,gyze,isign,mixup,sct,indx,nxe,nxh)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tfft = tfft + time
        !
        ! calculate magnetic field in fourier space with standard procedure:
        ! updates byze, wm
        call dtimer(dtime,itime,-1)
        call BBPOIS13(cue,byze,ffc,ci,wm,nx,nxeh,nxh)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tfield = tfield + time
        !
        ! transform magnetic force to real space with standard procedure:
        ! updates byze, gyze
        call dtimer(dtime,itime,-1)
        isign = 1
        call FFT1R2X(byze,gyze,isign,mixup,sct,indx,nxe,nxh)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tfft = tfft + time
        !
        ! add constant to magnetic field with standard procedure: updates byze
        call dtimer(dtime,itime,-1)
        call BADDEXT1(byze,omy,omz,nx,nxe)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tfield = tfield + time
        !
        ! copy guard cells with standard procedure: updates fxe, byze
        call dtimer(dtime,itime,-1)
        call DGUARD1L(fxe,nx,nxe)
        call CGUARD1L(byze,nx,nxe)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tguard = tguard + time
        !
        ! add longitudinal and old transverse electric fields with standard
        ! procedure: updates exyze
        call dtimer(dtime,itime,-1)
        call ADDVRFIELD13(exyze,cus,fxe,nxe)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tfield = tfield + time
        !
        ! deposit electron acceleration density and momentum flux with standard
        ! procedure: updates dcu, amu
        call dtimer(dtime,itime,-1)
        dcu = 0.0; amu = 0.0
        call GDJPOST1L(part,exyze,byze,dcu,amu,omx,qme,qbme,dt,idimp,np,  &
            &nxe)
        ! add old scaled electric field with standard procedure: updates dcu
        call ASCFGUARD1L(dcu,cus,q2m0,nx,nxe)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tdcjpost = tdcjpost + time
        !
        ! add guard cells with standard procedure: updates dcu, amu
        call dtimer(dtime,itime,-1)
        call ACGUARD1L(dcu,nx,nxe)
        call ACGUARD1L(amu,nx,nxe)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tguard = tguard + time
        !
        ! transform acceleration density and momentum flux to fourier space
        ! with standard procedure: updates dcu, amu, gyze
        call dtimer(dtime,itime,-1)
        isign = -1
        call FFT1R2X(dcu,gyze,isign,mixup,sct,indx,nxe,nxh)
        call FFT1R2X(amu,gyze,isign,mixup,sct,indx,nxe,nxh)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tfft = tfft + time
        !
        ! take transverse part of time derivative of current with standard
        ! procedure: updates dcu
        call dtimer(dtime,itime,-1)
        call ADCUPERP13(dcu,amu,nx,nxeh)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tfield = tfield + time
        !
        ! calculate transverse electric field with standard procedure:
        ! updates cus, wf
        call dtimer(dtime,itime,-1)
        isign = -1
        call EPOIS13(dcu,cus,isign,ffe,ax,affp,wpm,ci,wf,nx,nxeh,nxh)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tfield = tfield + time
        !
        ! transform transverse electric field to real space with standard
        ! procedure: updates cus, gyze
        call dtimer(dtime,itime,-1)
        isign = 1
        call FFT1R2X(cus,gyze,isign,mixup,sct,indx,nxe,nxh)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tfft = tfft + time
        !
        ! copy guard cells with standard procedure: updates cus
        call dtimer(dtime,itime,-1)
        call CGUARD1L(cus,nx,nxe)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tguard = tguard + time
        !
        ! add longitudinal and transverse electric fields with standard
        ! procedure: exyze = cus + fxe, updates exyze
        ! cus needs to be retained for next time step
        call dtimer(dtime,itime,-1)
        call ADDVRFIELD13(exyze,cus,fxe,nxe)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tfield = tfield + time
        !
        ! inner iteration loop
        do k = 1, ndc
            !
            ! deposit electron current and acceleration density and momentum flux
            ! with standard procedure: updates cue, dcu, amu
            call dtimer(dtime,itime,-1)
            cue = 0.0; dcu = 0.0; amu = 0.0
            call GDCJPOST1L(part,exyze,byze,cue,dcu,amu,omx,qme,qbme,dt,idimp,&
                &np,nxe)
            ! add scaled electric field with standard procedure: updates dcu
            call ASCFGUARD1L(dcu,cus,q2m0,nx,nxe)
            call dtimer(dtime,itime,1)
            time = real(dtime)
            tdcjpost = tdcjpost + time
            !
            ! add guard cells for current, acceleration density, and momentum flux
            ! with standard procedure: updates cue, dcu, amu
            call dtimer(dtime,itime,-1)
            call ACGUARD1L(cue,nx,nxe)
            call ACGUARD1L(dcu,nx,nxe)
            call ACGUARD1L(amu,nx,nxe)
            call dtimer(dtime,itime,1)
            time = real(dtime)
            tguard = tguard + time
            !
            ! transform current to fourier space with standard procedure:
            ! update cue, gyze
            call dtimer(dtime,itime,-1)
            isign = -1
            call FFT1R2X(cue,gyze,isign,mixup,sct,indx,nxe,nxh)
            call dtimer(dtime,itime,1)
            time = real(dtime)
            tfft = tfft + time
            !
            ! calculate magnetic field in fourier space with standard procedure:
            ! updates byze, wm
            call dtimer(dtime,itime,-1)
            call BBPOIS13(cue,byze,ffc,ci,wm,nx,nxeh,nxh)
            call dtimer(dtime,itime,1)
            time = real(dtime)
            tfield = tfield + time
            !
            ! transform magnetic force to real space with standard procedure:
            ! updates byze, gyze
            call dtimer(dtime,itime,-1)
            isign = 1
            call FFT1R2X(byze,gyze,isign,mixup,sct,indx,nxe,nxh)
            call dtimer(dtime,itime,1)
            time = real(dtime)
            tfft = tfft + time
            !
            ! add constant to magnetic field with standard procedure: updates bzye
            call dtimer(dtime,itime,-1)
            call BADDEXT1(byze,omy,omz,nx,nxe)
            call dtimer(dtime,itime,1)
            time = real(dtime)
            tfield = tfield + time
            !
            ! transform acceleration density and momentum flux to fourier space
            ! with standard procedure: updates dcu, amu, gyze
            call dtimer(dtime,itime,-1)
            isign = -1
            call FFT1R2X(dcu,gyze,isign,mixup,sct,indx,nxe,nxh)
            call FFT1R2X(amu,gyze,isign,mixup,sct,indx,nxe,nxh)
            call dtimer(dtime,itime,1)
            time = real(dtime)
            tfft = tfft + time
            !
            ! take transverse part of time derivative of current with standard
            ! procedure: updates dcu
            call dtimer(dtime,itime,-1)
            call ADCUPERP13(dcu,amu,nx,nxeh)
            call dtimer(dtime,itime,1)
            time = real(dtime)
            tfield = tfield + time
            !
            ! calculate transverse electric field with standard procedure:
            ! updates cus, wf
            call dtimer(dtime,itime,-1)
            isign = -1
            call EPOIS13(dcu,cus,isign,ffe,ax,affp,wpm,ci,wf,nx,nxeh,nxh)
            call dtimer(dtime,itime,1)
            time = real(dtime)
            tfield = tfield + time
            !
            ! transform transverse electric field to real space with standard
            ! procedure: updates cus, gyze
            call dtimer(dtime,itime,-1)
            isign = 1
            call FFT1R2X(cus,gyze,isign,mixup,sct,indx,nxe,nxh)
            call dtimer(dtime,itime,1)
            time = real(dtime)
            tfft = tfft + time
            !
            ! copy guard cells with standard procedure: updates byze, cus
            call dtimer(dtime,itime,-1)
            call CGUARD1L(byze,nx,nxe)
            call CGUARD1L(cus,nx,nxe)
            call dtimer(dtime,itime,1)
            time = real(dtime)
            tguard = tguard + time
            !
            ! add longitudinal and transverse electric fields with standard
            ! procedure: exyze = cus + fxyze, updates exyze
            ! cus needs to be retained for next time step
            call dtimer(dtime,itime,-1)
            call ADDVRFIELD13(exyze,cus,fxe,nxe)
            call dtimer(dtime,itime,1)
            time = real(dtime)
            tfield = tfield + time
            !
        enddo
        !
        !=================Edited by Jinxing Li======================
        ! perform potential diagnostic and unpack into array pott
        ! updates potc, pott, we
        call POTP1(qe,potc,ffc,we,nx,nxeh,nxh)
        call RDMODES1(potc,pott,nx,modesx,nxeh,modesxd)
        ! perform vector potential diagnostic and unpack into array vpott
        ! updates vpotc, vpott, wm
        call APOTP13(cue,vpotc,ffc,ci,wm,nx,nxeh,nxh)
        call RDVMODES1(vpotc,vpott,nx,modesx,2,nxeh,modesxd)
        !===========================================================
        ! push particles with standard procedure: updates part, wke
        wke = 0.0
        call dtimer(dtime,itime,-1)
        call GBPUSH13L(part,exyze,byze,omx,qbme,dt,dt,wke,idimp,np,nx,nxe,&
            &ipbc)
        call dtimer(dtime,itime,1)
        time = real(dtime)
        tpush = tpush + time
        !
        ! sort particles by cell for standard procedure
        if (sortime > 0) then
            if (mod(ntime,sortime)==0) then
                call dtimer(dtime,itime,-1)
                call DSORTP1XL(part,part2,npic,idimp,np,nx1)
                ! exchange pointers
                tpart => part
                part => part2
                part2 => tpart
                call dtimer(dtime,itime,1)
                time = real(dtime)
                tsort = tsort + time
            endif
        endif
        !
        if (ntime==0) then
            wt = we + wm
            write (*,*) 'Initial Total Field, Kinetic and Total Energies:'
            write (*,'(3e14.7)') wt, wke, wke + wt
            write (*,*) 'Initial Electrostatic, Transverse Electric and Mag&
                &netic Field Energies:'
            write (*,'(3e14.7)') we, wf, wm
        endif
        !====================Added by Jinxing Li====================
        !write to "para.txt" fiel
        if (ntime==0) then
            open(540,file=trim(adjustl(filename2)))
            write(540,"(9E12.4E2)"),real(np/Np_Interval),real(Np_Interval),real(nx),&
                &dt,real(T_Interval),real(nloop),omx,omy,1.0/ci
            write(540,"(9E12.4E2)"),vtx(1),vty(1),vx0(1),vtx(2),vty(2),vx0(2),0.0,0.0,0.0
            close(540)
        end if
        if (mod(ntime,T_Interval)==0) then
            call CURLF1(vpotc,Bw_k,nx,nxeh)
            call GRADF1(potc,Ex_k,nx,1,nxeh)
           ! do kk=1,nxe
           !     Ex_kr(kk)=exyze(1,kk)
           !     Ey_kr(kk)=exyze(2,kk)
           !     Ez_kr(kk)=exyze(3,kk)
           ! end do
           ! call FFT1RXX(Ex_kr, gxe,-1,mixup,sct,indx,nxe,nxh)
           ! call FFT1RXX(Ey_kr, gxe,-1,mixup,sct,indx,nxe,nxh)
           ! call FFT1RXX(Ez_kr, gxe,-1,mixup,sct,indx,nxe,nxh)

            !write to field.txt
            do kk=1,nxe
                write(100),exyze(1:3,kk),byze(1:2,kk)
            end do
            write(100),wke,we,wm,wf,we+wm

            !write Exk, Byk and Bzk into pot.txt; only kk=1 to modexsd=nx/4
            do kk=1,modesxd
                write(97),&
                    &real(Ex_k(kk)),  aimag(Ex_k(kk)),&
                    &real(Bw_k(1,kk)),aimag(Bw_k(1,kk)),&
                    &real(Bw_k(2,kk)),aimag(Bw_k(2,kk))
                ! real part, imag part
            end do
        endif
        !write 10 distribution to edist.dat
        if (mod(ntime,nloop_percent*10)==0) then
            call dist_statistics(part,np,nv,nth,vArr,thArr,vArr1,thArr1,fpsd,npsd,omx,omy)
            call dist_analytics(nkappa,ne,vtx,vty,vx0,kappa,nv-1,vArr1,fvpara,fvperp)
            do i=1,nv-1
                write(150,"(11E12.4E2)"),vArr1(i),&
                    &0.5*(fpsd(i,1)+fpsd(i,nth-1)),&
                    &0.5*(fpsd(i,2)+fpsd(i,nth-2)),&
                    &0.5*(fpsd(i,3)+fpsd(i,nth-3)),&
                    &0.5*(fpsd(i,4)+fpsd(i,nth-4)),&
                    &0.5*(fpsd(i,5)+fpsd(i,nth-5)),&
                    &0.5*(fpsd(i,6)+fpsd(i,nth-6)),&
                    &0.5*(fpsd(i,7)+fpsd(i,nth-7)),&
                    &fpsd(i,int(nth/2)),&
                    fvpara(i),fvperp(i)
            end do
        end if

        !===========================================================
    end do
    !ntime = ntime + 1
    !go to 500
    !2000 continue

    close(97)
    close(100)
    close(150)
    print *,"Successfully write to ",filename4
    print *,"============================="


    !
    ! * * * end main iteration loop * * *
    !
    write (*,*) 'ntime, ndc = ', ntime, ndc
    wt = we + wm
    write (*,*) 'Final Total Field, Kinetic and Total Energies:'
    write (*,'(3e14.7)') wt, wke, wke + wt
    write (*,*) 'Final Electrostatic, Transverse Electric and Magnetic&
        & Field Energies:'
    write (*,'(3e14.7)') we, wf, wm
    !
    write (*,*)
    write (*,*) 'deposit time = ', tdpost
    write (*,*) 'current deposit time = ', tdjpost
    write (*,*) 'current derivative deposit time = ', tdcjpost
    tdpost = tdpost + tdjpost + tdcjpost
    write (*,*) 'total deposit time = ', tdpost
    write (*,*) 'guard time = ', tguard
    write (*,*) 'solver time = ', tfield
    write (*,*) 'fft time = ', tfft
    write (*,*) 'push time = ', tpush
    write (*,*) 'sort time = ', tsort
    tfield = tfield + tguard + tfft
    write (*,*) 'total solver time = ', tfield
    time = tdpost + tpush + tsort
    write (*,*) 'total particle time = ', time
    wt = time + tfield
    write (*,*) 'total time = ', wt
    write (*,*)
    !
    wt = 1.0e+09/(real(nloop)*real(np))
    write (*,*) 'Push Time (nsec) = ', tpush*wt
    write (*,*) 'Deposit Time (nsec) = ', tdpost*wt
    write (*,*) 'Sort Time (nsec) = ', tsort*wt
    write (*,*) 'Total Particle Time (nsec) = ', time*wt
    !
    stop
end program
