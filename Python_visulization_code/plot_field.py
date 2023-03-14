#plot_field.py
import numpy as np
import matplotlib.pyplot as plt
import sys
dir1='/home/jli/results/pic/whistler_kappa92/'
#f=open(dir1+'para.dat','rb')
#M0=np.fromfile(file=f,dtype='f4').reshape((2,9))
#f.close()
M0=np.loadtxt(dir1+'para.txt',comments='#')
dt0=M0[0,3]
T_interval=M0[0,4]
dt=dt0*T_interval
nx=M0[0,2]
nxe=nx+2
nxe1=nxe+1
omx=M0[0,6]
omy=M0[0,7]
print(T_interval,dt,nxe1)
print(omy)
print("omy")

print("Np/Np_inter Np_int nx dt T_inter nloop omx omy omz")
print(M0[0,:])
print("vtx vty vtz vx0 vy0 vz0")
print(M0[1,:])
f=open(dir1+'field.dat','rb')
M=np.fromfile(file=f,dtype='f4').reshape((-1,5))
f.close()
#       0   1   2   3   4
#0      ex  ey  ez  by  bz
#...
#nxe+1  wke we  wm  wf  we+wm

# wke: particle kinetic energy
#  we: electrostatic energy
#  wf: magnetic field energy
#  wm: transverse electric field energy
#  wt: total energy

[len1,len2]=np.shape(M)
NT=int(len1/nxe1)
print(len1,nxe1,NT)

plt.figure(figsize=(12,9))

plt.subplot(2,2,1)
plt.plot(np.arange(NT)*dt,M[nxe::nxe1,1],'r',label='electrostatic')
plt.plot(np.arange(NT)*dt,M[nxe::nxe1,2],'g',label='magnetic')
plt.plot(np.arange(NT)*dt,M[nxe::nxe1,3],'b',label='E transverse')
plt.plot(np.arange(NT)*dt,M[nxe::nxe1,4],'k',label='total em')
plt.plot(np.arange(NT)*dt,M[nxe::nxe1,0]-M[nxe,0],'c',label='particle')
plt.xlabel('$T * \omega_{pe}$')
plt.ylabel('$Energy$')
plt.title('energy')
plt.legend()

plt.subplot(2,2,2)
Z=np.zeros([nxe,NT])
for k in range(0,NT):
	Z[:,k]=M[(nxe1*k):(nxe1*k+nxe),0]
#Z=np.transpose(Z)
plt.imshow(Z,cmap='jet',origin='lower',interpolation='none',extent=[0,NT*dt,0,nxe],aspect='auto')
plt.xlabel('$T * \omega_{pe}$')
plt.ylabel('$X /\lambda_{DE}$')
plt.title('Ex component')
plt.colorbar()

plt.subplot(2,2,3)
Z=np.zeros([nxe,NT])
for k in range(0,NT):
	Z[:,k]=M[(nxe1*k):(nxe1*k+nxe),3]-omy
#Z=np.transpose(Z)
bmax=np.max(Z)
bmin=np.min(Z)
if (bmax > -bmin):
	bmax=-bmin
else:
	bmin=-bmax

plt.imshow(Z,cmap='jet',origin='lower',interpolation='none',extent=[0,NT*dt,0,nxe],aspect='auto',vmin=bmin,vmax=bmax)
plt.xlabel('$T * \omega_{pe}$')
plt.ylabel('$X /\lambda_{DE}$')
plt.title('By component')
plt.colorbar()

plt.subplot(2,2,4)
Z2=np.zeros([nxe,NT])
for k in range(0,NT):
	Z2[:,k]=M[(nxe1*k):(nxe1*k+nxe),1]
#Z2=np.transpose(Z2)
print(np.shape(Z2))
plt.imshow(Z2,cmap='jet',origin='lower',interpolation='none',extent=[0,NT*dt,0,nxe],aspect='auto')
plt.xlabel('$T * \omega_{pe}$')
plt.ylabel('$X /\lambda_{DE}$')
plt.title('Ey component')
plt.colorbar()

plt.show()

