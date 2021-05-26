import numpy as np
import pylab
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors
from scipy import integrate 
import os
import sys
#for python 3.6

printtofile=True
showplots=False ## if set to True printtofile will fail

# define n, T, B profiles and plasma beta (dimensionless pressure)
te_power=-0.5
def beta(x): 
	#x in this function is heliospheric distance in units of AU
	beta=8*np.pi*(4*x**(-2))*(1e7*1.6022e-19*12.21*x**(te_power))/(10419.9e-5*(0.0233/x)**2*np.sqrt(1+(x)**2))**2 #define beta as 8 pi n T / B^2, and plugging in explicit expressions for n, T, B
	return beta

def dens(r):
	n=4*r**(-2)
	n=n*1e6 # convert to particles/m^3
	return n

def Bfield(r):
	r_45=1
	ro=0.0233
	Bo=10419.9e-5
	B=Bo*(ro/r)**2*np.sqrt(1+(r/r_45)**2) #B in units of Gauss 
	return B

Tempe=lambda r: 12.21*r**(te_power) #in eV
Tempi=lambda r: 12.21*r**(-0.7)

# take in system argument for radius and compute relevant quantities 
radius = 0.01*float(sys.argv[1]) #0.31 #1.0 # AU
angle = float(sys.argv[2])

print('radius in AU is ', radius)
densityatr= dens(radius) 
delta= Bfield(radius)**2/(0.137**2 * 1e-6*dens(radius))  # for 8nT magnetic field, delta=8.5e-8
Te=Tempe(radius)
Ti=Tempi(radius)
beta_e=beta(radius)
beta_i=(Ti/Te)*beta_e
print('elc beta at r is ',beta_e)
print('ion beta at r is ',beta_i)
print('~~~~~~~~~~~~~~~~~~~~~')
print('core elc temp at r in eV is ',Te)
print('core ion temp at r in eV is ',Ti)
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('density at r in particles/m^3 is ',densityatr)
print('B strength at r in G is ',Bfield(radius))
print('delta at r is ',delta)

# edit input file to change betai and delta values
linenum=30
file_name='input.dat'
fid = open(file_name).readlines()

d=fid[linenum-15].split()
fid[linenum-15]=fid[linenum-15].replace(d[2],str(angle))
d=fid[linenum-14].split()
fid[linenum-14]=fid[linenum-14].replace(d[2],str(delta))

d=fid[linenum].split()
fid[linenum]=fid[linenum].replace(d[2],str(beta_i))
d=fid[linenum+1].split()
fid[linenum+1]=fid[linenum+1].replace(d[2],str(beta_i))

with open((file_name),'w') as l:
	for j in fid:
		l.write(j)

#######################

Nspecies=2

#######################


npara=np.zeros(Nspecies,dtype='i4')
nperp=np.zeros(Nspecies,dtype='i4')

vparamin=np.zeros(Nspecies)
vparamax=np.zeros(Nspecies)

vperpmin=np.zeros(Nspecies)
vperpmax=np.zeros(Nspecies)

dens=np.zeros(Nspecies)
mu=np.zeros(Nspecies)
beta_para=np.zeros(Nspecies)
beta_perp=np.zeros(Nspecies)
vdrift=np.zeros(Nspecies)

########################

#species 1

npara[0]= 525
nperp[0]= 150

vparamin[0]=-125
vparamax[0]=325

vperpmin[0]=0.0
vperpmax[0]=125


#species 2

npara[1]=npara[0] # set velocity resolution equal to species 1 (elc strahl) so that they can be added together
nperp[1]=nperp[0]

mu[1]=9.1e-31 
beta_para[1]=10.0
beta_perp[1]=10.0

vparamin[1]=vparamin[0] # set velocity limits equal to species 1 (elc strahl) so that they can be added together
vparamax[1]=vparamax[0] 

vperpmin[1]=vperpmin[0] 
vperpmax[1]=vperpmax[0] 

#############################

limit=10.0**(-300)

def dist_bimax(vpar,vper, n,m,beta_par,beta_per,drift):
        bimax=np.exp(-n*(vpar-drift)**2/beta_par/m -n*vper**2/beta_per/m)* n**1.5 /(m**1.5 *np.pi**1.5 *beta_per*np.sqrt(beta_par))
        return bimax.ravel()

def dist_isotropic_max(vpar,vperp,m,drift,r=radius,temp=Te,delt=delta):
	#no_si = 4e6  # m^-3
	e = 1.6022e-19 #J/eV
	#T=lambda r: 12.21*r**(-0.5) #in eV -- this follows 1 AU temp from Horaites & Astfalk et al
	#Ao = (2*np.pi*T(r)*e/m)**(-3/2) 
	Ao = (2*np.pi*Te*e/m)**(-3/2)
	va=2.9866e8*np.sqrt(delta)
	maxwell=va**3*Ao*np.exp((-(vpar*va-drift*va)**2 - (vperp*va)**2 ) / (2*Te*e/m)) # norm by alf sp
	#maxwell=va**3*Ao*np.exp((-(vpar*va-drift*va)**2 - (vperp*va)**2 ) / (2*T(r)*e/m)) # norm by alf sp
	#maxwell=no_si*Ao*np.exp(-m*((vpar*va-drift*va)**2 + (vperp*va)**2) / (2*T(r)*e)) #norm by alf sp & density
	return maxwell.ravel()

def dist_strahl(vpar,vper,r=radius,delt=delta,diagnostics=False): #,n,m,beta_par,beta_per,drift):
	
	# convert (vpar,vperp) order pair to (v,mu) i.e. speed & cosine pitch angle
	v = np.sqrt(vpar**2 + vper**2)
	try:
		mu = vpar / v
	except ZeroDivisionError:
		mu = 1.0
	
	T = lambda r: 12.21*r**(te_power) # typical solar wind temperature scaling with Te=12.21eV at 1AU
	n = lambda r: 4e6*r**(-2) # typical solar wind density scaling with 4cm^-3 at 1AU 
	ro=0.0233 # 5 solar radii
	To = T(ro)
	no_si = n(ro)	
	
	e_phi_inf = 4*To # from Boldyrev, Egedal, Forest (2020); balance flux of ions and electrons at ro --> e*phi_inf/To ~ 4
	e = 1.6022e-19 # J/eV
	beta = 1.05 # this is (1+Zeff)/2
	Lambd = 20 # perhaps more appropriate for solar corona
	
	me_si = 9.1e-31 #kg
	eps_o = 8.854e-12 #F/m 
	
	Ao = no_si*(me_si/(2*np.pi*To*e))**(3/2) #normalize to density near corona
	lambdo = 4*np.pi*eps_o**2*(e*To)**2/(no_si*e**4*Lambd*beta) / 1.496e+11 #in units of AU
	e_phi_inf = 4*To
	Fo = (e_phi_inf/To)*np.exp(-e_phi_inf/To) 
	
	# constants here at 1 AU
	va=2.9866e8*np.sqrt(delta) # delta is squared ratio of cyclotron freq to plasma freq, used as parameter in LEOPARD input file
 	
	#define hyperbolic tangent cutoff function to avoid disconinutities
	cutoff=1.0
	C = lambda vpar: np.exp(cutoff*vpar)/(np.exp(cutoff*vpar)+1)
		
	# radial functions
	r_45 = 1 #approx location where magnetic field lines are at 45 deg angle, in units of AU
	B = lambda r: (ro/r)**2*np.sqrt(1+(r/r_45)**2)  # here B is B(r) normalized to Bo - here neglecting specification of Bo, just using Parker spiral scaling for field strength 

	KE = (me_si/2)*(v*va)**2/e # in eV
	deltaE = KE - T(r) #attempt to normalize by alfven speed

	R = lambda r: r*(1 - 2*(T(r)/deltaE) + 2*(T(r)/deltaE)**2*np.log(deltaE/T(r) + 1)) #eqn 11 in Boldyrev & Horaites 2019

	strahl =  C(vpar) * va**3*Ao*Fo * (lambdo/R(r)) * ((deltaE+e_phi_inf)/e_phi_inf)*(deltaE/To)*np.exp(-deltaE/To) * np.exp(-deltaE*KE*(1-mu**2)*lambdo/(To**2*R(r)*B(r))) #eqn 16/17 in Boldyrev 2019
		
	if diagnostics:	
		expfactor=-deltaE*KE*(1-mu**2)*lambdo/(To**2*R(r)*B(r))
		energy=deltaE*KE*(1-mu**2)/(To**2)
		dist=lambdo/R(r)
		mag=1/(B(r))
		print('---------------------------------')
		print('-------- velocities -------------')
		print('v parallel is ', vpar,' and v perp is ', vper)
		print('distance is ', r)
		print('-------- contants -----------')
		print('To in eV is ', To)
		print('no in cm^-3 is ', no_si*1e-6)
		print('Ao is ', Ao)
		print('Fo is ', Fo)
		print('mean free path in AU is ', lambdo)
		#print('cgs mean free path is ', lambdo_cgs)
		print('e_phi_inf in eV is ', e_phi_inf)
		print('alfven speed is ', va)
		print('kinetic energy alf normed is ', KE)
		print('kinetic energy non normed is ', me_si*v**2/(2*e))
		print('deltaE is ', deltaE)

		print('-------- exponetial factors ------')
		print('exponential factor is ', expfactor)
		print('energy dimensionless is ', energy)
		print('R function is', R(r))
		print('T/deltaE is ', T(r)/deltaE)
		print('T/deltEsqr*log(deltE/T +1) is ', 2*(T(r)/deltaE)**2*np.log(deltaE/T(r) + 1))
		print('dimensionless distance is ', dist)
		print('magnetic field factor is ', mag)
		print('total distribution is ', strahl)
		print('----------------------------------')	
	return strahl.ravel()

#print out densities and parallel strahl current for diagnostics
me_si=9.1e-31
mi_si=1.67e-27
va=2.9866e8*np.sqrt(delta)
print('alfven velocity is')
print(va*1e-3,' km/s')


dens[0]=1/densityatr # divide out density assigned in strahl function
print('density of strahl over requested grid')
print(integrate.dblquad(lambda y,x: 2*np.pi*x*dens[0]*dist_strahl(y,x), 0.0,vperpmax[0],lambda x:vparamin[0],lambda x:vparamax[0]))
#diagnostic extras 
#print('density of antisunward strahl over requested grid zero cutoff')
#print(integrate.dblquad(lambda y,x: 2*np.pi*x*dens[0]*dist_strahl(y,x), 0.0,vperpmax[0],lambda x:0.0,lambda x:vparamax[0]))
#print('density of sunward strahl over requested grid zero cutoff')
#print(integrate.dblquad(lambda y,x: 2*np.pi*x*dens[0]*dist_strahl(y,x), 0.0,vperpmax[0],lambda x:vparamin[0],lambda x:0.0))

print('parallel current of strahl')
print(integrate.dblquad(lambda y,x: 2*np.pi*x*y*dens[0]*dist_strahl(y,x), 0.0, vperpmax[0],lambda x:vparamin[0],lambda x:vparamax[0]))

# set density and parallel current of core correspondingly
dens[1]=1-integrate.dblquad(lambda y,x: 2*np.pi*x*dens[0]*dist_strahl(y,x), 0.0,vperpmax[0],lambda x:vparamin[0],lambda x:vparamax[0])[0]
vdrift[1]= -integrate.dblquad(lambda y,x: 2*np.pi*x*y*dens[0]*dist_strahl(y,x), 0.0, vperpmax[0],lambda x:vparamin[0],lambda x:vparamax[0])[0] / dens[1]
print('drift of core is',vdrift[1])
print('density of core over requested grid')
print(integrate.dblquad(lambda y,x: 2*np.pi*x*dens[1]*dist_isotropic_max(y,x,m=me_si,drift=vdrift[1]), 0.0,vperpmax[1],lambda x:vparamin[1],lambda x:vparamax[1]))
#parallel current should have same magnitude but opposite sign of strahl parallel current
print('parallel current of core over requested grid')
print(integrate.dblquad(lambda y,x: 2*np.pi*x*y*dens[1]*dist_isotropic_max(y,x,m=me_si,drift=vdrift[1]), 0.0,vperpmax[1],lambda x:vparamin[1],lambda x:vparamax[1]))

for ispecies in range(0,Nspecies):

	#file_name='distribution'+np.str(ispecies+1)+'.dat'

	vpara = np.linspace(vparamin[ispecies],vparamax[ispecies], npara[ispecies])
	vperp = np.linspace(vperpmin[ispecies],vperpmax[ispecies], nperp[ispecies])

	vpara2,vperp2=np.meshgrid(vpara,vperp)
	
	if ispecies==0:
		data=dist_strahl(*(vpara2, vperp2))
	else:
		data=dist_isotropic_max(*(vpara2,vperp2),m=mu[ispecies],drift=vdrift[ispecies],r=radius,delt=delta)
	data=data.reshape(nperp[ispecies],npara[ispecies])

	data_new=np.zeros((nperp[ispecies],npara[ispecies]))

	for i in range(0,npara[ispecies]):
		for j in range(0,nperp[ispecies]):

			if(data[j,i]>limit):
				data_new[j,i]=data[j,i]*dens[ispecies]
			else:
				data_new[j,i]=0.0

	if ispecies==0:
		data0=data_new
	elif ispecies==1:
		data1=data_new
#save total electron dist to file
ispecies=0  
vpara = np.linspace(vparamin[0],vparamax[0], npara[0])
vperp = np.linspace(vperpmin[0],vperpmax[0], nperp[0])
file_name='distribution1.dat'
dat_fin = open(file_name, 'w')
for i in range(0,npara[ispecies]):
	for j in range(0,nperp[ispecies]):
		dat_fin.write(str(vpara[i]))
		dat_fin.write(" ")
		dat_fin.write(str(vperp[j]))
		dat_fin.write(" ")
		dat_fin.write(str(data0[j,i]+data1[j,i]))
		dat_fin.write("\n")
if printtofile:
	#os.system('mv distribution1.dat {}/distribution/'.format(sys.argv[1]))
	os.system('mv distribution1.dat distribution/')

### extra diagnostics: plot electron distributions
vpara = np.linspace(vparamin[0],vparamax[0], npara[0])
vperp = np.linspace(vperpmin[0],vperpmax[0], nperp[0])
vpara2,vperp2=np.meshgrid(vpara,vperp)

#2d plot

#######################################################
# plot strahl only

plt.figure(2)
plt.pcolormesh(vpara,vperp,data0,cmap='inferno',norm=colors.LogNorm(vmin=1e-60,vmax=(data0).max()))
plt.colorbar()
plt.xlabel('$v_\Vert /\ v_A$')
plt.ylabel('$v_\perp /\ v_A$')
plt.title('electron strahl distribution')


plt.figure(3,figsize=(6,3.7))
plt.pcolormesh(vpara,vperp,data0+data1,cmap='inferno',norm=colors.LogNorm(vmin=1e-24,vmax=(data0+data1).max()))
plt.colorbar()
plt.xlabel('$v_\Vert /\ v_A$')
plt.ylabel('$v_\perp /\ v_A$')
plt.title('electron core + strahl distribution')
plt.axis([-160,350,0,160])
#plt.title('$f_s v_A^3$')
#######################################################

#3d plot
fig = plt.figure(0)
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(vpara2, vperp2, data0+data1,cmap='inferno')


plt.figure(4)
#np.size(vpara)
plt.plot(vpara,data0[0,:]+data1[0,:])
plt.xlabel('$v_\Vert /\ v_A$')
plt.ylabel('$f(v_\Vert, v_\perp = 0) v_A^3 $')

plt.figure(5)
plt.pcolormesh(vpara,vperp,data0,cmap='inferno')
plt.colorbar()
plt.xlabel('$v_\Vert /\ v_A$')
plt.ylabel('$v_\perp /\ v_A$')
plt.title('electron strahl distribution')

#display plots
if showplots:
	plt.show() 
