import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import math
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['ps.fonttype']=42
hfont={'fontname':'Helvetica'}
font= {'size':7}
matplotlib.rc('font',**font)

dir1='/home/jli/results/pic/whistler_kappa50/'
M=np.loadtxt(dir1+'edist.txt',comments='#')
M0=np.loadtxt(dir1+'para.txt',comments='#')
omx=M0[0,6]
omy=M0[0,7]
cc=20.0#M0[0,8]
fpefce=1.0/np.sqrt(omx**2+omy**2)
print("fpe/fce = ", fpefce)

#  vArr  fvpara   fvperp  nvpara  nvperp fvpara_ana  fvperp_ana
[len1,len2]=np.shape(M)
nv=int(np.around(len1/10))
print(len1,len2,nv)

plt.figure(figsize=(6,4.5))
fig,axes=plt.subplots(nrows=2,ncols=2)

plt.subplot(2,2,1)
ind1=np.arange(nv)
gammaArr=1.0/np.sqrt(1.0-(M[ind1,0]/cc)**2)
ekArr=(gammaArr-1.0)*9.1e-31*9.0e16/1.6e-19*1e-3
mksz=1
lw=0.5
plt.plot(ekArr,M[ind1,8],'.-r',linewidth=lw,markersize=mksz,label='PA=84-96 deg')#fvpara_statistics
plt.plot(ekArr,M[ind1,6],'.-g',linewidth=lw,markersize=mksz,label='PA=60-72 deg')#fvpara_statistics
plt.plot(ekArr,M[ind1,4],'.-c',linewidth=lw,markersize=mksz,label='PA=36-48 deg')#fvperp_statistics
plt.plot(ekArr,M[ind1,1],'.-k',linewidth=lw,markersize=mksz,label='PA=0-12 deg')#fvpara_statistics
#plt.plot([1.0/fpefce/2.0,1.0/fpefce/2.0],[1e-9,1e4],':c')
plt.xlabel('Energy (keV)')
plt.ylabel('PSD')
plt.yscale('log')
plt.xscale('log')
leg=plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.ylim(1e-8,1e1)
plt.xlim(1e-2,3e2)

#------------------------- Subplot 2 ------------------
plt.subplot(2,2,2)
ind2=np.arange(nv)+nv*6
gammaArr=1.0/np.sqrt(1.0-(M[ind2,0]/cc)**2)
ekArr=(gammaArr-1.0)*9.1e-31*9.0e16/1.6e-19*1e-3

plt.plot(ekArr,M[ind2,8],'.-r',linewidth=lw,markersize=mksz,label='PA=84-96 deg')#fvpara_statistics
plt.plot(ekArr,M[ind2,6],'.-g',linewidth=lw,markersize=mksz,label='PA=60-72 deg')#fvpara_statistics
plt.plot(ekArr,M[ind2,4],'.-c',linewidth=lw,markersize=mksz,label='PA=36-48 deg')#fvperp_statistics
plt.plot(ekArr,M[ind2,1],'.-k',linewidth=lw,markersize=mksz,label='PA=0-12 deg')#fvpara_statistics
#plt.plot([1.0/fpefce/2.0,1.0/fpefce/2.0],[1e-9,1e4],':c')
print(cc/fpefce/2.0)
print(cc)
plt.xlabel('Energy (keV)')
plt.ylabel('PSD')
plt.yscale('log')
plt.xscale('log')
leg=plt.legend()
leg.get_frame().set_linewidth(0.0)

plt.ylim(1e-8,1e1)
plt.xlim(1e-2,3e2)



pi=3.1415927
paras=np.loadtxt(dir1+'para.txt',comments='#')
f=open(dir1+'pot.dat','rb')
pot=np.fromfile(file=f,dtype='f4').reshape((-1,6))
f.close()
# Exr,Exi,Byr,Byi,Bzr,Bzi

nx  =int(round(paras[0,2]))
dnp =int(round(paras[0,1])) #Interval of particle when output
npx =int(round(paras[0,0]))
T_interval=int(round(paras[0,4]))#T Interval when output
dt0=paras[0,3]*T_interval #dt in this program

nloop_tot=int(round(paras[0,5])) # Total time series
omx =paras[0,6]
omy =paras[0,7]
omz =paras[0,8]

nxe=nx+2
nmod=int(round(nx/4))
[len1,len2]=np.shape(pot)
nt=round(len1/nmod)
print("txt.y  nt")
print(len1,nt)

#----------------------------Subplot 3 -------------------
pot1=pot[round(len1/100*5):round(len1/100*13),:]
[len3,len4]=np.shape(pot1)
nt=round(len3/nmod)

Fs=1.0/dt0;
kArr=np.arange(nmod-1)-0.5
wArr=2.0*pi*Fs*np.arange(round(nt/2))/nt/np.sqrt(omx**2+omy**2) 
dw=2.0*pi*Fs/nt/omx;
nf=np.size(wArr);
wk=np.zeros((nf,nmod));
print(nloop_tot,T_interval,nt)
print(nmod,nf)


for k in range(0,nmod):
	SampleArr=pot1[k:((nt-1)*nmod+k+1):nmod,5 ]
	LenSample=np.size(SampleArr)
	SampleArr0 = SampleArr[0]+np.arange(LenSample)/LenSample*(SampleArr[-1]-SampleArr[0])
	SampleArr  =SampleArr-SampleArr0
	Hann=np.sin(pi*np.arange(LenSample)/(LenSample-1))**2;
	SampleArr  =SampleArr*Hann;
	Y1=np.fft.fft(SampleArr)
	Y2=abs(Y1)
	Y3=2.0*Y2[0:round(nt/2) ]
	Y3[0]=Y3[0]*0.5
	wk[:,k]=np.log10(Y3)


plt.subplot(2,2,3)
plt.imshow(wk,cmap='jet',origin='lower',interpolation='none',extent=[kArr[0],kArr[-1],wArr[0],wArr[-1]],aspect='auto',vmin=-2.2,vmax=-0.2)#,norm=LogNorm())
plt.plot([kArr[0],kArr[-1]],[0.5,0.5],'--r',linewidth=lw);
#plt.plot([kArr[0],kArr[-1]],[1.0,1.0],'--r',linewidth=lw);

plt.xlabel('k mode')
plt.ylabel('$\omega\: / \:\omega_{ce}$')
#plt.colorbar()
plt.ylim(0.05,1.05)
plt.xlim(0.0,22.0)

#----------------------------- Subplot 4 --------------------------
pot1=pot[round(len1/100*60):round(len1/100*68),:]
[len3,len4]=np.shape(pot1)
nt=round(len3/nmod)

Fs=1.0/dt0;
kArr=np.arange(nmod-1)-0.5
wArr=2.0*pi*Fs*np.arange(round(nt/2))/nt/np.sqrt(omx**2+omy**2) 
dw=2.0*pi*Fs/nt/omx;
nf=np.size(wArr);
wk=np.zeros((nf,nmod));
print(nloop_tot,T_interval,nt)
print(nmod,nf)


for k in range(0,nmod):
	SampleArr=pot1[k:((nt-1)*nmod+k+1):nmod,5 ]
	LenSample=np.size(SampleArr)
	SampleArr0 = SampleArr[0]+np.arange(LenSample)/LenSample*(SampleArr[-1]-SampleArr[0])
	SampleArr  =SampleArr-SampleArr0
	Hann=np.sin(pi*np.arange(LenSample)/(LenSample-1))**2;
	SampleArr  =SampleArr*Hann;
	Y1=np.fft.fft(SampleArr)
	Y2=abs(Y1)
	Y3=2.0*Y2[0:round(nt/2) ]
	Y3[0]=Y3[0]*0.5
	wk[:,k]=np.log10(Y3)


plt.subplot(2,2,4)
plt.imshow(wk,cmap='jet',origin='lower',interpolation='none',extent=[kArr[0],kArr[-1],wArr[0],wArr[-1]],aspect='auto',vmin=-2.2,vmax=-0.2)#,norm=LogNorm())
plt.plot([kArr[0],kArr[-1]],[0.5,0.5],'--r',linewidth=lw);
#plt.plot([kArr[0],kArr[-1]],[1.0,1.0],'--r',linewidth=lw);


plt.xlabel('k mode')
plt.ylabel('$\omega\: / \:\omega_{ce}$')
#plt.colorbar()
plt.ylim(0.05,1.05)
plt.xlim(0.0,22.0)
xxx=[1,5,10,15,20]
xx1=range(0,4)
fig.tight_layout()

plt.savefig('/home/jli/results/CaseStudy/chorus_gap/aa.pdf',format='pdf')
