import numpy as np
import matplotlib.pyplot as plt

path='' #### path of your folder

# box parameters
L = 60.0  # length box [0, L]
N = 100  # Number of pixels
d = L / N
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)

X, Y = np.meshgrid(x, y)

Nt = 300  # time
dt = 0.3



hbar = 1.0
m = 1.0   

D = hbar*1j/(2*m)


def Psi_vector(Psi_mat):
    i=0
    j=0
    NL,Nl=np.shape(Psi_mat)
    psi = np.zeros((NL**2), dtype=complex)
    for k in range(NL**2):

        psi[k]=Psi_mat[i,j]
        i+=1
        if i==NL:
            i=0
            j+=1

    return psi


def Psi_matrice(Psi_vec):
    i=0
    j=0
    l=len(Psi_vec)
    psi=np.zeros((int(np.sqrt(l)),int(np.sqrt(l))), dtype=complex)
    for k in range(l):
        psi[i,j]=Psi_vec[k]

        i+=1

        if i==int(np.sqrt(l)):
            i=0
            j+=1

    return psi


# Potential
V = np.zeros((N,N), complex) 

# latice parameters
spacing_Y = 5
spacing_X = 5
number_Y = 5
number_X = int(L/spacing_X)
R_atome = 0.5
for i in range(number_X-1):
    X_p = (i+1)*spacing_X*N/L
    for j in range(number_Y):
        Y_p = (N/2) + j*spacing_Y*N/L
       
        for u in range(N):
            for v in range(N):
                if np.sqrt((u - X_p)**2 + (v - Y_p)**2) <= R_atome:
                    V[u, v] = 100

Vvector = Psi_vector(V)

#latice visualization
plt.imshow(np.abs(V), extent=(0, L, 0, L), origin='lower', aspect='auto',cmap='hot')
plt.show()


# initial wave function parameter
k0x = -2
k0y = 0
delta_kx = 0.2
delta_ky = 0.2
x0 = L/4
y0 = L/2

# initial  wave function
psi0 = np.sqrt(delta_kx * delta_ky / np.pi) * \
    np.exp(-(((delta_kx * (X - x0))**2) / 2) + 1j * k0x * (X - x0)) * \
    np.exp(-(((delta_ky * (Y - y0))**2) / 2) + 1j * k0y * (Y - y0))

# Display of the initial probability density
density0 = np.real(np.conj(psi0)*psi0)
fig, ax = plt.subplots()
plt.imshow(density0, extent=(0, L, 0, L), origin='lower', aspect='auto',cmap='hot')
plt.colorbar()
plt.ylim(0,L)
plt.xlim(0,L)
plt.savefig(path+"/"+"0000"+".png")
plt.show()
plt.close('all')



print("processing...")
# Cranck-Nicolson
alpha = D * dt / (2 * d**2)

A =  np.diag(np.ones(N**2)) + np.diag(4 * alpha * np.ones(N**2))
A += np.diag(-alpha * np.ones(N**2-1), 1) + np.diag(-alpha * np.ones(N**2-1), -1)
A += np.diag(-alpha * np.ones(N**2-N), N) + np.diag(-alpha * np.ones(N**2-N), -N)
A += np.diag(1j*dt*Vvector/2)

B =  np.diag(np.ones(N**2)) + np.diag(-4 * alpha * np.ones(N**2)) 
B += np.diag(alpha * np.ones(N**2-1), 1) + np.diag(alpha * np.ones(N**2-1), -1)
B += np.diag(alpha * np.ones(N**2-N), N) + np.diag(alpha * np.ones(N**2-N), -N) 
B += np.diag(-1j*dt*Vvector/2)

AinvB = np.linalg.solve(A,B)


Psi = Psi_vector(psi0)

file_n=1

for n in range(Nt):   
    Psi= np.dot(AinvB, Psi)

    Psi_trace = Psi_matrice(Psi)

    Density_proba = np.real(np.conj(Psi_trace)*Psi_trace)

    fig, ax = plt.subplots()
    plt.imshow(Density_proba, extent=(0, L, 0, L), origin='lower', aspect='auto',cmap='hot')
    plt.colorbar()
    plt.text(0.05, 1.01*L, 't = '+str(np.around((n+1)*dt,3)))
    plt.ylim(0,L)
    plt.xlim(0,L) 
    plt.savefig(path+"/"+("0"*(4-len(str(file_n))))+str(file_n)+".png")
    print(str(file_n)+"/"+str(Nt))
    file_n += 1
    plt.close('all')


import os
from moviepy.editor import ImageSequenceClip

print(os.getcwd())
image_folder = path

os.chdir(image_folder)
 
images = [img for img in os.listdir(image_folder)
        if  img.endswith(".jpg") or
            img.endswith(".jpeg") or
            img.endswith("png")]
     
print(images) 
  
clip = ImageSequenceClip(images, fps = 40)

clip.ipython_display(width = 360)

