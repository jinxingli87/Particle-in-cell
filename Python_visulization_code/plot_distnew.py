import numpy as np
import matplotlib.pyplot as plt
import sys
import math
dir1='/home/jli/results/pic/whistler_kappa45/'
#f=open(dir1+'edist.dat','rb')
#M=np.fromfile(file=f,dtype='f4').reshape((-1,7))
#f.close()
M=np.loadtxt(dir1+'edist.txt',comments='#')
M0=np.loadtxt(dir1+'para.txt',comments='#')
omx=M0[0,6]
omy=M0[0,7]
cc=M0[0,8]
fpefce=1.0/np.sqrt(omx**2+omy**2)
print("fpe/fce = ", fpefce)

#  vArr  fvpara   fvperp  nvpara  nvperp fvpara_ana  fvperp_ana
[len1,len2]=np.shape(M)
#print(M[:,0])
nv=int(np.around(len1/10))
print(len1,len2,nv)

plt.figure(figsize=(12,9))
plt.subplot(1,2,1)
ind1=np.arange(nv)+nv*0
plt.plot(M[ind1,0],M[ind1,8],'.-r',label='fvperp')
plt.plot(M[ind1,0],M[ind1,6],'.-g',label='fvpara')
plt.plot(M[ind1,0],M[ind1,3],'.-b',label='fvpara')
plt.plot(M[ind1,0],M[ind1,1],'.-k',label='fvpara')
plt.plot(M[ind1,0],M[ind1,9],':k',label='fvpara_ana')#fvpara_analytic
plt.plot(M[ind1,0],M[ind1,10],':r',label='fvperp_ana')#fvperp_analytic
plt.plot([cc/fpefce/2.0,cc/fpefce/2.0],[1e-9,1e4],':c')
plt.xlabel('Vecolity / vt')
plt.ylabel('PSD')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.ylim(1e-9,1e4)
plt.xlim(1e-1,20.0)


plt.subplot(1,2,2)
ind2=np.arange(nv)+round(nv*2)
#ind2=np.arange(nv*9,nv*10)
plt.plot(M[ind1,0],M[ind1,1],':k',label='fvpara')
plt.plot(M[ind1,0],M[ind1,8],':r',label='fvperp')

plt.plot(M[ind2,0],M[ind2,8],'.-r',label='fvperp')
plt.plot(M[ind2,0],M[ind2,6],'.-g',label='fvpara')
plt.plot(M[ind2,0],M[ind2,3],'.-b',label='fvpara')
plt.plot(M[ind2,0],M[ind2,1],'.-k',label='fvpara')
plt.plot([cc/fpefce/2.0,cc/fpefce/2.0],[1e-9,1e4],':c')
print(cc/fpefce/2.0)
print(cc)
plt.xlabel('Vecolity / vt')
plt.ylabel('PSD')
plt.yscale('log')
plt.xscale('log')
plt.legend()

plt.ylim(1e-9,1e4)
plt.xlim(1e-1,20.0)
plt.show()
