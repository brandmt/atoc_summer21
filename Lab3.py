import numpy as np
import matplotlib.pyplot as plt
import miepython as mp

ri = 0.3 #microns
s = 2 #sigma
ni = (1e4)*1e-12 #µm^-3
r = np.arange(0.01,10.01,0.01) #microns

def dndr(r):
    dn = ni*np.exp(-((np.log(r/ri))**2/(2*(np.log(s))**2)))
    dr = np.sqrt(2*np.pi)*r*np.log(s)
    return dn/dr

N = dndr(r)*0.01 #Number Distribution
wl = 550 #nm
wl = wl/1000 #µm
xa = 2*np.pi*r/wl #Size Parameter

#Refractive Indices
aerosols = ['Water Soluble','Dustlike','Soot','Oceanic']
rfr = np.array([1.53,1.53,1.75,1.381])
rfi = np.array([6e-3,8e-3,4.4e-1,4.26e-9])

for i in range(4):
    ma  = complex(rfr[i],-rfi[i])
    qext, qsca, qback, g = mp.mie(ma,xa)
    
    Bext = np.sum(qext*np.pi*r**2*N)
    Bsca = np.sum(qsca*np.pi*r**2*N)
    Babs = Bext - Bsca
    scat_alb = Bsca/Bext
      
    print(aerosols[i])
    print('Extinction =',round(Bext*1e9,2))
    print('Scattering =',round(Bsca*1e9,2))
    print('Absorption =',round(Babs*1e9,10))
    print('Scattering Albedo =',round(scat_alb,10),'\n')
 
#%% Problem 3
wlnm = np.array([400,600,700,900,1060,1500])
wls = wlnm/1000
Bext = np.empty(6)
plt.figure(3,dpi=200)
for i in range(4):
    ma  = complex(rfr[i],-rfi[i])
    for j in range(6):
        xa = 2*np.pi*r/wls[j]
        wave = str(wlnm[j])
        qext, qsca, qback, g = mp.mie(ma,xa)
        Bext[j] = np.sum(N*qext*np.pi*r**2)
    plt.plot(wlnm,Bext*1e9,label=aerosols[i])
    plt.scatter(wlnm,Bext*1e9,s=10)
    plt.grid(alpha=0.2)
    plt.ylabel('Extinction Coefficient $(km^{-1}$)')
    plt.xlabel('Wavelength (nm)')
    plt.legend()
 
#%% Problem 5
aerosols = ['Water Soluble','Oceanic']
r = np.arange(0.1,10.01,0.01)
rfr = np.array([1.53,1.381])
rfi = np.array([6e-3,4.26e-9])
theta = np.array([60,90])

plt.figure(4,figsize=(12,10),dpi=200)
for i in range(2): #Two Angles
    mu = np.array([np.cos(np.deg2rad(theta[i]))])
    plt.subplot(2, 1, i+1)
    plt.title('\u03B8='+str(theta[i])+u'\N{DEGREE SIGN}',fontweight='bold')
    for j in range(2): #Two Aerosols
        ma  = complex(rfr[j],-rfi[j])
        sca_int = np.empty(len(r))
        for k in range(len(r)): #Each Radii
            xa = np.pi*r[k]/wl
            sca_int[k] = mp.i_unpolarized(ma, xa, mu)
        plt.plot(r,sca_int,label=aerosols[j])
    plt.grid(alpha=0.2)
    plt.ylabel('Scattering Intensity $(sr^{-1}$)')
    plt.xlabel('Radius (µm)')
    plt.legend()
    # plt.xscale('log')
