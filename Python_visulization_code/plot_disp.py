# plot_disp.py
# plot the wave power spectral density on (w,k)
import numpy as np
import matplotlib.pyplot as plt
#plot_disp.py
import sys
from matplotlib.colors import LogNorm

pi=3.1415927
dir1='/home/jli/results/pic/whistler_kappa45/'
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

pot1=pot[round(len1/100*25):round(len1/100*30),:]
pot1=pot[round(len1/100*3):round(len1/100*8),:]
pot1=pot[round(len1/100*0):round(len1/100*5),:]
pot1=pot[round(len1/100*0):round(len1/100*10),:]
pot1=pot[round(len1/100*40):round(len1/100*45),:]
[len1,len2]=np.shape(pot1)
nt=round(len1/nmod)


Fs=1.0/dt0;
kArr=np.arange(nmod-1)
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
	Y2=abs(Y1)**2
	Y3=2.0*Y2[1:round(nt/2)+1 ]
	#Y3[0]=Y3[0]*0.5
	wk[:,k]=np.log10(Y3)

plt.figure(figsize=(12,9))

plt.imshow(wk,cmap='jet',origin='lower',interpolation='none',extent=[kArr[0],kArr[-1],wArr[0],wArr[-1]],aspect='auto',vmin=-4.5,vmax=-0.2)#,norm=LogNorm())
plt.plot([kArr[0],kArr[-1]],[0.5,0.5],'--r');
plt.plot([kArr[0],kArr[-1]],[1.0,1.0],'--r');

plt.xlabel('k mode')
plt.ylabel('Frequenc / wce')
plt.colorbar()
plt.ylim(0.0,1.2)
plt.xlim(0.0,20.0)


plt.show()
sys.exit()
'''
[X,Y]=ndgrid(kArr,wArr);

hf=figure(1);
set(hf,'Position',[100,100,800,600]);
%set(hf,'PaperPosition',[100,100,450,700])
%set(hf, 'PaperPositionMode', 'auto')   % Use screen size
%set(hf, 'PaperUnits', 'points'); % or centimeters
colormap('jet');
tmp=log(wk);
tmp0=min(min(tmp));
tmp1=max(max(tmp));

h2=pcolor(X,Y,tmp);
caxis([tmp1-8.0,tmp1]);
colorbar;


hold on;
plot([0,max(kArr)],[0.5,0.5],'--w');
plot([0,max(kArr)],[1.0,1.0],'--w');
%set(gca,'Yscal','log');
set(h2,'EdgeColor','none');
%set(gca,'fontsize',FS);
set(gca,'Layer','top');
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
FS=14;
set(gca,'TickDir','out','fontsize',FS);
xlabel('k mode','fontsize',FS)
ylabel('Frequency (\omega/\omega_{ce})','fontsize',FS);
caxis([-3,2.5]);

end

'''
