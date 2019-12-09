import numpy as np
from math import pi, sqrt, sin, cos, radians, tan, asin, exp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from random import randrange
import tqdm
import random
from numpy.linalg import norm
from scipy import special
from astropy import constants as const
from scipy.special import wofz


#Parametros para el modelo

dnub = 100 #diametro de las nubes en pc
radis = 30000 #radio del disco en pc
hdis = 10000 #altura del disco pc

xygrill =  int((radis/dnub))
nh = int((hdis/dnub))

volnub = (4/3)*pi*(dnub/2)**3

D = np.load('r_mage_kpc.npy')

def probabilidad(r):
    return(9*10**(-8)*r**2-0.006*r+100)

'''radplot = np.linspace(1, 30000, 1000)
radplot1 = np.linspace(1,30,100)
probplot = []
velplot = []

for i in range(len(radplot)):
    probplot.append(probabilidad(radplot[i]))

plt.plot(radplot, probplot)
plt.show()'''

#Velocidades
inc = radians(45) #radianes -----> obs
alpha = radians(14) #radianes ----> obs????
hv = 10 #kpc  -----> parametro libre
vR = 200 #km/s ----> parametro libre

velobs = np.load('vel1.npy')



def vel(r, caca):
    x0 = abs(caca) * cos(alpha) * 1000
    print(x0)
    y = sqrt((r**2)- (x0)**2) #pc
    vrot = sqrt(1+(sin(alpha)**2)*tan(inc)**2) * velobs[5] /(sin(inc)*cos(alpha))
    vang = (vrot * sin(inc) * exp(-abs(int(y))/(hv*tan(inc)))/sqrt(1+(y/(D[5]*1000)**2)))
    vr = vR * sin(inc) * y /sqrt((y**2)+(((D[5]*1000)**2)))
    return(vang+vr)


'''for i in range(len(radplot1)):
    velplot.append(vel(radplot1[i]))

plt.plot(radplot1, velplot)
plt.show()'''

def numnubestrue(impactparam):
    D1 = impactparam*1000
    xc = D1 * cos(alpha)
    y0 = D1 * sin(alpha)/cos(inc)
    z = np.linspace(-hdis/2, hdis/2, nh+1)
    print(z)
    n = 0
    velos = []
    radios = []
    for i in range(len(z)):
        zc = z[i]
        if zc<=0:
            yc = (-1) * (((zc+1000)/(2*tan(inc)))-y0)
        else:
            yc = (((zc+1000)/(2*tan(inc)))+y0)

        rc = sqrt(xc**2 + yc**2)
        veli = vel(rc, impactparam)
        prob = probabilidad(rc)

        selec = random.uniform(0, 100)
        if selec<prob:
            n = n+1
            velos.append(veli)
            radios.append(rc)
        else:
            pass
    return(n,velos,radios)

#calcular tau

N = 10**(12)
redshift = 0.73379
b = 10

lam0, f, gamma, mass = [2796.3521, 0.6155, 2.68e8, 24.305]
c  = const.c.to('cm/s').value
sigma0 = 0.0263

def Tau(lam,vel):

    lamc = ((vel/const.c.to('km/s').value)+1)*((lam0))
    nu = c/(lam*1e-8)
    nu0 = c/(lamc*1e-8)

    dnu = nu - (nu0/(1+redshift))
    dnud = (b*100000)*nu/c

    x = dnu/dnud
    y = gamma/(4*pi*dnud)
    z = x + 1j*y
    v = np.real(wofz(z)/(sqrt(pi)*dnud))
    tau = N * sigma0 * f  *v
    vt = np.exp(-tau)
    obs_spec_wave = lam / (1+redshift)
    dv = (const.c.to('km/s').value * ((lam / (lamc *(1 + redshift)))-1))

    velt = []
    for i in range(len(dv)):
        velt.append(dv[i]+vel)

    return(tau,dv)


#Suma tau
fig1, axs1 = plt.subplots(1, 4)
fig2, axs2 = plt.subplots(1,4)
fig3, axs3 = plt.subplots(1,3)

Ns = []

sampling = 0.005 #en armstrong
lamini  = 4835
lamfin = 4900
numlam = (lamfin-lamini)/sampling

wejemplo = np.linspace(lamini,lamfin,numlam)


for i in range(len(D)):
    Ns.append(numnubestrue(D[i]))

Fluxes = []
for j in range(len(Ns)):
    velejemplo = Ns[j][1]
    if len(velejemplo)== 0:
        pass
    else:
       #wejemplo = np.linspace(4835,4850,1000)

       tausum = []
       velsum = []
       for i in range(len(velejemplo)):
           tau = Tau(wejemplo, velejemplo[i])
           vele = tau[1]
        #print(vele)
           if j<4:
               axs1[j].step(vele,np.exp(-tau[0]),'b')
               axs1[j].set_ylim(-0.05,1.05)
               axs1[j].set_title('%s, %s' %(j, Ns[j][0]))
           if 3<j<8:
                axs2[(j-4)].step(vele,np.exp(-tau[0]), 'b')
                axs2[(j-4)].set_ylim(-0.05,1.05)
                axs2[j-4].set_title('%s, %s' %(j, Ns[j][0]))
           if 7<j<12:
               axs3[(j-8)].step(vele,np.exp(-tau[0]),'b')
               axs3[(j-8)].set_ylim(-0.05,1.05)
               axs3[j-8].set_title('%s, %s' %(j, Ns[j][0]))
           else:
                pass

           tausum.append(tau[0])

       tausum = np.asarray(tausum)
    #print(tausum)
       tauejemplo=np.sum(tausum,0)
       Flux = np.exp(-tauejemplo)
       Fluxes.append(Flux)
       #Fluxes = Fluxes[:-1]
       if j<4:
           axs1[j].step(vele, Flux, 'r')
           axs1[j].set_ylim(-0.05,1.05)
       if 3<j<8:
           axs2[(j-4)].step(vele,Flux, 'r')
           axs2[(j-4)].set_ylim(-0.05,1.05)
       if 7<j<12:
           axs3[(j-8)].step(vele,Flux, 'r')
           axs3[(j-8)].set_ylim(-0.05,1.05)
       else:
           pass

plt.show()


vejemplo = (const.c.to('km/s').value * ((wejemplo / (lam0 *(1 + redshift)))-1))

plt.figure()
poster1 = []
for i in range(len(Ns[1][1])):
    taup = Tau(wejemplo, Ns[1][1][i])
    plt.plot(vejemplo,np.exp(-taup[0]))
    poster1.append(taup[0])

tautp=np.sum(poster1,0)
fluxplot = np.exp(-tautp)
np.savetxt('velteo.dat', vejemplo)
np.savetxt('fluxteo.dat', Fluxes)
plt.plot(vejemplo, np.exp(-tautp), 'r')
plt.xlabel('vel[km/s]')
plt.ylabel('Norm. Flux')
plt.xlim(-10,200)
plt.ylim(-0.05, 1.05)
plt.show()



  #Filtro Gauss

from astropy.convolution import convolve, Gaussian1DKernel

fw = 3 #pixel

gauss_kernel = Gaussian1DKernel(3)

gausfluxt = []

fig4, axs4 = plt.subplots(1, 4)
fig5, axs5 = plt.subplots(1,4)
fig6, axs6 = plt.subplots(1,3)


for i in range(len(Fluxes)):
    gausflux = convolve(Fluxes[i], gauss_kernel)
    gausfluxt.append(gausflux)

for i in range(len(gausfluxt)):
           if i<4:
               axs4[i].step(vele,gausfluxt[i],'b')
               axs4[i].set_ylim(-0.05,1.05)
               axs4[i].set_title('%s, %s' %(i, Ns[i][0]))
           if 3<i<8:
                axs5[(i-4)].step(vele,gausfluxt[i], 'b')
                axs5[(i-4)].set_ylim(-0.05,1.05)
                axs5[i-4].set_title('%s, %s' %(i, Ns[i][0]))
           if 7<i<12:
               axs6[(i-8)].step(vele,gausfluxt[i],'b')
               axs6[(i-8)].set_ylim(-0.05,1.05)
               axs6[i-8].set_title('%s, %s' %(i, Ns[i][0]))
           else:
                pass


    #plt.xlim(4000,6000)
    #plt.ylim(0.1)
plt.show()



velocidadesplot = []
radiosplot = []

for i in range(len(Ns)):
    N = Ns[i]
    for j in range(len(N[1])):
        velocidadesplot.append(N[1][j])
        radiosplot.append(N[2][j])

plt.plot(radiosplot,velocidadesplot)
plt.show()
