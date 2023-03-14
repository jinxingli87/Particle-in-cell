# plot_spec.py
# plot the wave power spectral density on (w,t)
import numpy as np
import matplotlib.pyplot as plt
#plot_disp.py
import sys
from matplotlib.colors import LogNorm

pi=3.1415927
dir1='/home/jli/results/pic/whistler_kappa50/'
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
nt0=round(len1/nmod)
print("txt.y  nt")
print(len1,nt0)

#nmod=256
#nt0=50000
nshift=500*nmod#round(len1/MT)# nshift=2% T_total
nbox=nshift*2
nf=round(nbox/nmod/2)#this is because each episode is 4% of the total samples
nt_spec=round((len1-nbox)/nshift)+1
Fs=1.0/dt0;
ww=np.zeros((nf,nt_spec))
ww1=ww
print("nf, nk, nt_spec")
print([nf,nmod,nt_spec])

for n in range(0,nt_spec):
	pot1=pot[nshift*n:(nshift*n+nbox),:]
	[len1,len2]=np.shape(pot1)
	nt=round(len1/nmod)
	wArr=2.0*pi*Fs*np.arange(round(nt/2))/nt/np.sqrt(omx**2+omy**2) 
	dw=2.0*pi*Fs/nt/omx;
	nf=np.size(wArr);

	kArr=np.arange(nmod-1)
	wk=np.zeros((nf,nmod));

	for k in range(0,nmod):
		SampleArr=pot1[k:((nt-1)*nmod+k+1):nmod,5 ]
		LenSample=np.size(SampleArr)
		SampleArr0 = SampleArr[0]+np.arange(LenSample)/LenSample*(SampleArr[-1]-SampleArr[0])
		SampleArr  =SampleArr-SampleArr0
		Hann=np.sin(pi*np.arange(LenSample)/(LenSample-1))**2;
		SampleArr  =SampleArr*Hann;
		Y1=np.fft.fft(SampleArr)
		Y2=abs(Y1)**2
		Y3=2.0*Y2[1:round(nt/2)+1 ]
		#Y3[0]=Y3[0]*0.5
		wk[:,k]=Y3

	for m in range(0,nf):
		ww[m,n]=sum(wk[m,2:20])

plt.figure(figsize=(12,9))

plt.imshow(np.log10(ww),cmap='jet',origin='lower',interpolation='none',extent=[0,nt0*dt0,wArr[0],wArr[-1]],aspect='auto',vmin=-4.0,vmax=0.0)#,norm=LogNorm())

plt.xlabel('Time *wpe')
plt.ylabel('Frequenc / wce')
plt.colorbar()
plt.ylim(0.0,1.2)
plt.xlim(0.0,nt0*dt0)


plt.show()
sys.exit()
