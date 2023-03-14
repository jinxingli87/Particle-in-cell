subroutine edist(part,nk,nx,np,ne,vtx,vty,vx0,kappa,log_vmin,log_vmax,nvel,omx,omy)
    implicit none
    real x,y,z
    integer i,j,k
    real gaussian_dist,kappa_dist
    real,dimension(:,:),pointer:: dvdv,dist2d
    real,dimension(:),pointer:: vArr,vArr1,dist1d
    real,dimension(4) :: sum_part
    real ne,vtx,vty,vx0,kappa,ne1
    integer nvel1,nvel,nk
    dimension ne(nk),vtx(nk),vty(nk),vx0(nk),kappa(nk),ne1(nk)
    real sum1,sum2,sum3,tmp,tmp1,tmp2
    integer np,nx
    real part
    dimension part(4,np)
    double precision randum1
    real log_dv
    real log_vmin,log_vmax,omx,omy,bk,bl,sinth,costh
    integer i_tmp

    nvel1=nvel-1
    log_dv=(log_vmax-log_vmin)/nvel1
    allocate(vArr(nvel),vArr1(nvel1))
    allocate(dvdv(nvel1,nvel1),dist2d(nvel1,nvel1),dist1d(nvel1*nvel1))
    dist2d=0.0
    dist1d=0.0

    sum1=0.0
    ! Calculate vArr and dvArr
    !vArr    1    2    3    4 ...   nV
    !vArr1      1   2    3     nV-1

    ! Generate velocity array
    do i=1,nvel
        vArr(i)=10**(log_dv*real(i-1)+log_vmin)
    end do
    do i=1,nvel1
        vArr1(i)=0.5*(vArr(i)+vArr(i+1))
    end do
    vArr(1)=0.0

    do j=1,nk
        ne1(j)=real(ne(j))/real(sum(ne))
    end do

    ! Case of Kappa distribution
    do i=1,nvel1
        do j=1,nvel1
            dvdv(i,j)=(vArr(i+1)-vArr(i)) * (vArr(j+1)-vArr(j))
            do k=1,nk
                dist2d(i,j)=dist2d(i,j)+ne1(k)*0.25*&
                &(kappa_dist(kappa(k),vtx(k),vty(k),vx0(k),vArr(i),vArr(j))&
                &+kappa_dist(kappa(k),vtx(k),vty(k),vx0(k),vArr(i+1),vArr(j))&
                &+kappa_dist(kappa(k),vtx(k),vty(k),vx0(k),vArr(i),vArr(j+1))&
                &+kappa_dist(kappa(k),vtx(k),vty(k),vx0(k),vArr(i+1),vArr(j+1)))
            end do
        end do
    end do
    !v_para
    !  |
    !  |
    ! i|(i-1)*Nv+1 (i-1)*Nv+2 .. (i-1)*Nv+j ...
    !  |
    !  |
    ! 2|      Nv+1       Nv+2 ..       Nv+j ..
    ! 1|         1          2  ...        j ...
    !  __________________________________________
    !        j=  1          2  ...        j            v_perp

    do i=1,nvel1
        do j=1,nvel1
            sum1=sum1+dvdv(i,j)*dist2d(i,j)*2.0*3.14159*(vArr(j)+vArr(j+1))
            dist2d(i,j)=sum1 
            ! now the dist2d is the integrated possibility, and is
            ! monocromaticlly increasing from 0 towards 1
        end do
    end do
    print *,"----------------2"
    print *,"sum= ", sum1 !sum1 should be close to 1 but less than 1.
    
    do i=1,nvel1
        do j=1,nvel1
            dist1d(j+(i-1)*nvel1)=dist2d(i,j)/sum1
        end do
    end do

    ! now the dist1d is the integrated possibility
    
!       vperp
! vpara  1,1    1,2    1,3 ...
!        2,1    2,2    2,3 ...
!        3,1    3,2    3,3 ...
! (1,1)->1;  (2,1)->2; (3,1)->3

    tmp1=(10.0**log_dv-1.0)*0.45 !c ratio of variation in velocity
    do i=1,np
        part(1,i)=real(nx)/real(np)*(real(i)-0.5)
        tmp=real(randum1())
        j=1
        do while (tmp > dist1d(j))
            j=j+1
        end do
        part(2,i)=vArr1(int(j/nvel1)+1)*(1.0+tmp1*(2.0*real(randum1())-1.0))&
            &*sign(1.0,real(randum1())-0.5)  !c vpara
        tmp=real(randum1()) ! vperp phase
        tmp2=2.0*real(randum1())-1.0 ! random number from -1 to 1
        i_tmp=mod(j,nvel1)
        if (i_tmp==0) then
            i_tmp=nvel1
        end if
        part(3,i)=vArr1(i_tmp)*(1.0+tmp1*tmp2)*cos(tmp*6.28318)
        part(4,i)=vArr1(i_tmp)*(1.0+tmp1*tmp2)*sin(tmp*6.28318)
    end do

    !ensure <v_para> and <v_perp> are zero if no beams
    sum_part=sum(part,2)/real(np)
    tmp=0.0
    do i=1,nk
        tmp=tmp+vx0(i)*ne(i)!tmp is the averaged beam velocity
    end do

    do i=1,np
        part(2,i)=part(2,i)-sum_part(2)+tmp
        part(3,i)=part(3,i)-sum_part(3)
        part(4,i)=part(4,i)-sum_part(4)
    end do

    !  v_perp   |y      /
    !\          |      / v_para,B
    !  \        |     /
    !    \      |    /
    !      \    |   /
    !        \  |  /  theta
    !          \|______________________ x,k
    !   
    !
    !coordinate transformation
    costh=omx/sqrt(omx**2+omy**2)
    sinth=omy/sqrt(omx**2+omy**2)
    do i=1,np
        bk= part(2,i)*costh-part(3,i)*sinth
        bl= part(2,i)*sinth+part(3,i)*costh
        part(2,i)=bk
        part(3,i)=bl
    end do

    !ensure <vx>,<vy> and <vz> are all zero except the beam
    sum_part=sum(part,2)/real(np)
    do i=1,np
        part(2,i)=part(2,i)-sum_part(2)+tmp*costh
        part(3,i)=part(3,i)-sum_part(3)+tmp*sinth
        part(4,i)=part(4,i)-sum_part(4)
    end do


end


subroutine dist_analytics(nk,ne,vtx,vty,vx0,kappa,nv,vArr,fvpara,fvperp)
    implicit none
    integer nk,nv,i,j
    real ne1(nk)
    integer ne(nk)
    real vtx(nk),vty(nk),vx0(nk),kappa(nk),vArr(nv),fvperp(nv),fvpara(nv)
    real kappa_dist

    fvpara=0.0
    fvperp=0.0

    ! Generate f(vpara) and f(vperp)
    do j=1,nk
        ne1(j)=real(ne(j))/real(sum(ne))
    end do

    do i=1,nv
        do j=1,nk
            fvpara(i)=fvpara(i)+ne1(j)*&
                &kappa_dist(kappa(j),vtx(j),vty(j),vx0(j),vArr(i),0.0)
            fvperp(i)=fvperp(i)+ne1(j)*&
                &kappa_dist(kappa(j),vtx(j),vty(j),vx0(j),0.0,vArr(i))
        end do
    end do
end


subroutine dist_statistics(part,np,nv,nth,vArr,thArr,vArr1,thArr1,fpsd,npsd,omx,omy)
    implicit none
    integer np,nv,nth,i,j,k
    real part,vArr,thArr,fpsd,npsd
    dimension part(4,np),vArr(nv),thArr(nth)
    dimension fpsd(nv-1,nth-1),npsd(nv-1,nth-1)
    real vArr1,dvArr,thArr1,dthArr
    dimension vArr1(nv-1),dvArr(nv-1),thArr1(nth-1),dthArr(nth-1)
    real v_tot,vperp,vpara,pa
    real omx,omy,sinth,costh
    ! Calculate vArr1 and dvArr
    !vArr    1    2    3    4 ...   nV
    !vArr1      1   2    3     nV-1
    do i=1,nv-1
        vArr1(i)=sqrt(vArr(i)*vArr(i+1))
        dvArr(i)=vArr(i+1)-vArr(i)
    end do

    ! Calculate thArr and dthArr
    do i=1,nth-1
        thArr1=0.5*(thArr(i)+thArr(i+1))
        dthArr=thArr(i+1)-thArr(i)
    end do

    !  v_perp   |y      /
    !\          |      / v_para,B
    !  \        |     /
    !    \      |    /
    !      \    |   /
    !        \  |  /  theta
    !          \|______________________ x,k
    !   
    !
    !coordinate transformation
    !coordinate transformation
    costh=omx/sqrt(omx**2+omy**2)
    sinth=omy/sqrt(omx**2+omy**2)

    npsd=0.0
    do i=1,np
        vpara = part(2,i)*costh+part(3,i)*sinth
        vperp =-part(2,i)*sinth+part(3,i)*costh
        vperp =sqrt(part(4,i)**2+vperp**2)!v_perp
        v_tot =sqrt(part(2,i)**2+part(3,i)**2+part(4,i)**2)!v_tot
        pa=atan(vperp/vpara)
        if (vpara < 0) then
            pa=pa+3.1415926
        end if

        j=1
        do while(v_tot > vArr(j+1) .AND. j<nv)
            j=j+1
        end do
        k=1
        do while(pa > thArr(k+1) .AND. k<nth)
            k=k+1
        end do
        npsd(j,k)=npsd(j,k)+1
    end do

    do j=1,nv-1 
        do k=1,nth-1
            fpsd(j,k)=npsd(j,k)/(6.2831853/3.0*(vArr(j+1)**3-vArr(j)**3)*&
                &(cos(thArr(k))-cos(thArr(k+1))))/np
        end do 
    end do

    
end




function kappa_dist(kappa,vtx,vty,vx0,vx,vy)
    implicit none
    real kappa,vtx,vty,vx0,vx,vy
    real gaussian_dist,kappa_dist
    real tmp1,tmp2,tmp3,tmp4
    real PI
    if (kappa ==0) then
        kappa_dist=gaussian_dist(vtx,vty,vx0,vx,vy)
    else
        PI=3.14159
        tmp1=2.0**(2*kappa-1)*(kappa-0.5) * gamma(kappa)**2
        tmp2=PI**2*sqrt(kappa)*gamma(2.0*kappa)*vtx*vty**2
        tmp3=1.0 + vy**2/vty**2/kappa + (vx-vx0)**2/vtx**2/kappa
        kappa_dist=tmp1/tmp2*tmp3**(-kappa-1.0)
    end if
    return
end function

function gaussian_dist(vtx,vty,vx0,vx,vy)
    implicit none
    real*4 vtx,vty,vx0,vx,vy,gaussian_dist
    real*4 PI
    PI=3.14159
    gaussian_dist=1.0/(PI**1.5*vty**2*vtx)*exp(-vy**2/vty**2-(vx-vx0)**2/vtx**2)
    return
end function




function randum1()
    ! this is a version of the random number generator dprandom due to
    ! c. bingham and the yale computer center, producing numbers
    ! in the interval (0,1).  written for the sun by viktor k. decyk, ucla
    implicit none
    integer isc,i1,r1,r2
    double precision randum1,h1l,h1u,r0,r3,asc,bsc
    save r1,r2,h1l,h1u
    data r1,r2 /1271199957,1013501921/
    data h1l,h1u /65533.0d0,32767.0d0/
    isc = 65536
    asc = dble(isc)
    bsc = asc*asc
    i1 = r1 - (r1/isc)*isc
    r3 = h1l*dble(r1) + asc*h1u*dble(i1)
    i1 = r3/bsc
    r3 = r3 - dble(i1)*bsc
    bsc = 0.5d0*bsc
    i1 = r2/isc
    isc = r2 - i1*isc
    r0 = h1l*dble(r2) + asc*h1u*dble(isc)
    asc = 1.0d0/bsc
    isc = r0*asc
    r2 = r0 - dble(isc)*bsc
    r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
    isc = r3*asc
    r1 = r3 - dble(isc)*bsc
    randum1 = (dble(r1) + dble(r2)*asc)*asc
    return
end function
