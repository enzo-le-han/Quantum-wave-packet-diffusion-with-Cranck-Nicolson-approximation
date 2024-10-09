import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from matplotlib.patches import Rectangle

path='' #### path of your folder

# box parameters
L = 60.0  # length box [0, L]
N = 90  # Number of pixels
d = L / N
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)

X, Y = np.meshgrid(x, y)

Nt = 300  # time
dt = 0.3 


hbar= 1.0
m = 1.0   

D = hbar*1j/(2*m)




def wall():
    color="w"
    if slit == 1:
        top_wall = Rectangle((j0*d,0),     Lf, i3*d,      color=color, alpha=0.3)
        mid_wall = Rectangle((j0*d,i2*d),  Lf,(i1-i2)*d,  color=color, alpha=0.3)
        bottom_wall    = Rectangle((j0*d,i0*d),  Lf, i3*d,      color=color, alpha=0.3)
        ax.add_patch(bottom_wall)
        ax.add_patch(mid_wall)
        ax.add_patch(top_wall)

    color="b"
    if border ==1:
        border_left = Rectangle((0,0),     d, N*d,      facecolor='gray', hatch='x')
        border_bottom    = Rectangle((0,0), N*d, d,      facecolor='gray',hatch='x')
        border_right    = Rectangle(((N-1)*d,0), d, N*d,      facecolor='gray',hatch='x')
        border_top    = Rectangle((0,(N-1)*d), N*d, d,      facecolor='gray',hatch='x')
        ax.add_patch(border_left)
        ax.add_patch(border_bottom)
        ax.add_patch(border_right)
        ax.add_patch(border_top)
  

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




# Double slit parameters
Lf = 2 * d  # Width of the walls of the double slit.
s = 6  # Length of the wall between the slits.
a = 2  # Length of the slits.

# Horizontal
j0 = int(1/(2*d)*(L-Lf))  # left border.
j1 = int(1/(2*d)*(L+Lf))  # right border.

# Vertical
i0 = int(1 / (2 * d) * (L + s) + a / d)   # Lower boundary of the lower slit.
i1 = int(1 / (2 * d) * (L + s))           # Upper boundary of the lower slit.
i2 = int(1 / (2 * d) * (L - s))           # Lower boundary of the upper slit.
i3 = int(1 / (2 * d) * (L - s) - a / d)   # Upper boundary of the upper slit.




# Potential
V = np.zeros((N,N), complex) 

# slit
slit = 1 # 0:no   1:yes
if slit == 1:
    potential_in_slit = 100
    V[0:i3, j0:j1] = potential_in_slit
    V[i2:i1,j0:j1] = potential_in_slit
    V[i0:,  j0:j1] = potential_in_slit

# border
border = 1 # 0:no   1:yes
if border == 1:
    potential_in_borderure = 50000
    V[:,0] = V[:,N-1] = V[0,:] = V[N-1,:] = potential_in_borderure

Vvector = Psi_vector(V)




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
area = simps(simps(density0, x), y)
fig, ax = plt.subplots()
plt.imshow(density0, extent=(0, L, 0, L), origin='lower', aspect='auto',cmap='hot')
plt.colorbar()
wall()
plt.text(0.05, 1.01*L, 't = '+str(0))
plt.text(0.05, 1.05*L, 'Normalization condition= '+str(np.around(area,decimals=9)))
plt.ylim(0,L)
plt.xlim(0,L)
plt.savefig(path+"/"+"0000"+".png")
plt.show()
plt.close('all')




# Cranck-Nicolson
print("processing ...")
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
    Normalization_condition =np.around(simps(simps(Density_proba, x), y),decimals=9)

    fig, ax = plt.subplots()
    plt.imshow(Density_proba, extent=(0, L, 0, L), origin='lower', aspect='auto',cmap='hot')
    plt.colorbar()
    wall()
    plt.text(0.05, 1.01*L, 't = '+str(np.around((n+1)*dt,3)))
    plt.text(0.05, 1.05*L, 'Normalization condition = '+str(Normalization_condition))
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
