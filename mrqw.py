"""
This is a simualtion of a quantum walk on a one-dimensional lattice of length n.
parameters controlling the walk are thetam etc. x is the position of the particle
p is the momenttum, and m is the qubit  memory of the lattice that records
a particle visit. 
The probability of particle being at position x at time step t is given by the 
evolution of the state vector v=kronecker(x,p,m) through the operator E. E=A.R.M, 
where A is the advection (translation of x by the momentum p), R is the 
action of memory m on the momentum p, and M is the interaction between memory
the particle to record the visit. All of these are unitary operators, defined in 
the paper "History Dependent Quantum Random Walks as Quantum Lattice Gas Automaa",
by Shakeel, Meyer, and Love.
"""
import numpy as np
from scipy.sparse import csc_matrix
from scipy.linalg import norm
import numpy as np

import  scipy.io
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt



n=11

thetam=0#np.pi/2.0
theta0=np.pi/4.0
theta1=np.pi/4.0
theta2=np.pi/4.0
theta3=np.pi/4.0

rlist=[] #row indices of the matrix. This is where the particle is
clist=[] # column indices of the matrix. This is where the particle will be.
vlist=[] # value of at the row, col index, i.e, the matrix coefficient.
cnt=0

# Creating the unitary matrix that multiplies the state of the particle at each 
# step
for xc in range(n): #
    for pc in range(2):
        for mc in range(2**n):
                bc=("{0:0"+str(n)+"b}").format(mc)[::-1]
                bcx=bc[xc:xc+1]
                bcxf=''.join('1' if t == '0' else '0' for t in bcx)
                bcf=bc[:xc]+bcxf+bc[xc+1:]
                mcf=int(bcf,2)
                rlist.append(xc*2*(2**n)+pc*(2**n)+mc)
                clist.append(xc*2*(2**n)+pc*(2**n)+mc)
                mrval= np.cos(thetam) #if mr == 0 else np.sin(np.pi/4.)*(1.0j)
                vlist.append(mrval)
                rlist.append(xc*2*(2**n)+pc*(2**n)+mcf)
                clist.append(xc*2*(2**n)+pc*(2**n)+mc)
                mrfval= np.sin(thetam)*(1.0j) #if mr == 0 else np.cos(np.pi/4.)
                vlist.append(mrfval) 
row  = np.array(rlist)
col  = np.array(clist)
data = np.array(vlist)
M=csc_matrix((data, (row, col)), shape=(n*2*(2**n), n*2*(2**n)), dtype=np.complex)


rlist=[]
clist=[]
vlist=[]
cnt=0
#thetab=np.pi/4.0
for xc in range(n):
    for pc in range(2):
        for mc in range(2**n):
                bc=("{0:0"+str(n)+"b}").format(mc)[::-1]
                bcx=bc[xc-1]+bc[np.fmod(xc+1,n)]
                mcx=int(bcx,2)
                if mcx==0:
                    thetab=theta0
                elif mcx==1:
                    thetab=theta1
                elif mcx==2:
                    thetab=theta2
                else:
                    thetab=theta3
                pcf= 0 if pc == 1 else 1
                rlist.append(xc*2*(2**n)+pc*(2**n)+mc)
                clist.append(xc*2*(2**n)+pc*(2**n)+mc)
                pcval= np.cos(np.pi/4.) if bcx == '0' else np.cos(thetab)
                vlist.append(pcval)
                rlist.append(xc*2*(2**n)+pcf*(2**n)+mc)
                clist.append(xc*2*(2**n)+pc*(2**n)+mc)
                pcfval=  np.sin(np.pi/4.)*(1.0j) if bcx == '0' else np.sin(thetab)*(1.0j)
                vlist.append(pcfval) 
row  = np.array(rlist)
col  = np.array(clist)
data = np.array(vlist)
R=csc_matrix((data, (row, col)), shape=(n*2*(2**n), n*2*(2**n)), dtype=np.complex)



rlist=[]
clist=[]
vlist=[]
cnt=0
for xc in range(n):
    for pc in range(2):
        for mc in range(2**n):
                #bc=("{0:0"+str(n)+"b}").format(mc)
                #bcx=bc[xc:xc+1]
                pci= -1 if pc == 1 else 1
                rlist.append((np.fmod(xc+pci+n,n))*2*(2**n)+pc*(2**n)+mc)
                clist.append(xc*2*(2**n)+pc*(2**n)+mc)
                vlist.append(1) 
row  = np.array(rlist)
col  = np.array(clist)
data = np.array(vlist)
A=csc_matrix((data, (row, col)), shape=(n*2*(2**n), n*2*(2**n)), dtype=np.complex)

E=A.dot(R.dot(M))


x=np.zeros(n, dtype=np.complex)
x[(n-1)/2]=1


p=np.zeros(2,dtype=np.complex)
p[0]=np.cos(np.pi/4.0)
p[1]=np.sin(np.pi/4.0)

m0=np.zeros(2,dtype=np.complex)
m0[0]=1
m1=np.zeros(2,dtype=np.complex)
m1[1]=1

m=m0
for k in range((n-1)/2-1):
    m=np.kron(m0,m)
m=np.kron(m0,m)
for k in range((n-1)/2):
    m=np.kron(m0,m)

    
v=csc_matrix(np.kron(x,np.kron(p,m)))
v=v.T   


Tsteps=(n-1)/2
pl=np.zeros((Tsteps+1,n))
for l in range(n):
    pl[0,l]=norm(v.todense()[l*(2**(n+1)):(l+1)*(2**(n+1)),0])
for k in range(Tsteps):
    #v=M.dot(v)
    #v=R.dot(v)
    #v=A.dot(v)
    v=E.dot(v)
    for l in range(n):
        pl[k+1,l]=np.around((norm(v.todense()[l*(2**(n+1)):(l+1)*(2**(n+1)),0]))**2, decimals=3) 


print pl   
    
plt.close("all")
fig = plt.figure()
fig.clf()
ax = fig.add_subplot(111, projection='3d')


# plotting the particle wave funtion as a function of time steps.
cnt=0
for c, z in zip(['r', 'g', 'c', 'b','m','y'][::-1], np.arange((n-1)/2+1)):
    xs = np.arange(-5,6)
    ys = pl[cnt,:]
    cnt=cnt+1  


    cs = [c] * len(xs)
    ax.bar(xs, ys, zs=z, zdir='y', color=cs, alpha=0.8)
ax.set_xlabel('position: x')
ax.set_ylabel('Time step: t')
ax.set_zlabel('Probability')
ax.set_title('Probability distribution of particle position', y=1.03)
ax.set_xlim3d(-((n-1)/2+1),(n-1)/2+1)
ax.set_ylim3d(0,(n-1)/2+1)
ax.set_zlim3d(0,1)
ax.set_ylim(ax.get_ylim()[::-1])
ax.view_init(azim=-50)
plt.show()

#plt.savefig('/Users/Asif/Documents/work/latex/ndqw2.eps', format='eps', dpi=100)

