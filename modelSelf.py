from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import rcParams 
import scipy.special as ss
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib import gridspec
import warnings
warnings.filterwarnings("ignore")
rcParams["font.size"] = 13
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
import emcee, corner
from numpy import matrix
import VBBinaryLensingLibrary as vb
VBB=vb.VBBinaryLensing()
VBB.Tol=1.0e-7; ### accuracy
VBB.SetLDprofile(VBB.LDlinear);
VBB.LoadESPLTable("./ESPL.tbl"); 
################################################################################
###    Constant           ###
G= 6.67430*pow(10.0,-11.0)
AU=1.495978707*pow(10.0,11)
Msun=1.989*pow(10.0,30.0)
Rsun =6.9634*pow(10.0,8.0)
KPC= 3.0857*pow(10.0,19.0)## meter 
velocity=299792458.0##m/s
KBol= 1.380649*pow(10.0,-23.0)
Hplank=6.62607015*pow(10.0,-34.0)
const= np.sqrt(4.0*G*Msun)/velocity
const2=float(Hplank*velocity*np.power(10.0,9.0)/KBol )




tp=0.0##in period
Mstar=2.08
Tstar=7542.0
fblend=1.0
Period=11.393423143561519#[days]
Dl=float(1.0/0.978796586561136)*KPC## meter
Rstar=3.8; 
################################################################################
Nw=110
throu=np.zeros((Nw , 2))
throu=np.loadtxt("./files/TESS_Throught.txt")
Nd=int(2603) 
tim=   np.zeros((Nd));  
Flux=  np.zeros((Nd));  
errf=  np.zeros((Nd));
par=np.zeros((Nd,3))
par=np.loadtxt("./files/LastData.dat")
tim,Flux,errf=par[:,0], par[:,1], par[:,2] 
tim=tim/Period#[-0.5, +0.5]
thre=float(0.001)
Nb=1000

################################################################################

NDIM=4+2
Nwalkers=16
Nstep= 300000
p0=np.zeros((Nwalkers,NDIM))
Nm=int(3346)
phi=   np.zeros((Nm));  
ksi=   np.zeros((Nm));  
x1=    np.zeros((Nm));  
x0=    np.zeros((Nm));
y1=    np.zeros((Nm));  
y0=    np.zeros((Nm));
z1=    np.zeros((Nm));  
Vx=    np.zeros((Nm));
RE=    np.zeros((Nm));
phase= np.zeros((Nm));
Astar= np.zeros((Nm));
term=  np.zeros((Nm));  
term0= np.zeros((Nm)); 
rho=   np.zeros((Nm));
u=     np.zeros((Nm));
dis=   np.zeros((Nm));
ampl=  np.zeros((Nm));
Delf1= np.zeros((Nm));
Delf2= np.zeros((Nm));
tmod=  np.zeros((Nm));
Tmod=  np.zeros((Nm));
for i in range(Nm):  
    tmod[i]=float(-0.53645364537465134+i/Nm/1.0034536451264523564) #[-0.5 , 0.5]
    phi[i] =(tmod[i]-tp)*2.0*np.pi 
################################################################################

def Fplank(wave):
    wave=wave*pow(10.0,-9.0)#[m]  
    con= Hplank*velocity/(KBol*Tstar*wave)
    return(2.0*Hplank*velocity*velocity*pow(wave,-5.0)/(np.exp(con)-1.0)  )  

def Doppler(vx):
    for m in range(Nm): 
        DF0=0.0
        Delf1[m]=0.0
        for s in range(Nw-1):
            wave= throu[s,0]             
            dw=float(throu[s+1,0]-throu[s,0])*pow(10.0,-9)
            waven= wave+wave*vx[m]/velocity
            Delf1[m]+= float(Fplank(waven) - Fplank(wave))*throu[s,1]*dw       
            DF0+=Fplank(wave) * throu[s,1] * dw
        Delf1[m]=float(Delf1[m]/DF0)     
    return(Delf1)

################################################################################

def Thirdlaw(mblak):
    return(np.power(G*Msun*(Mstar+mblak)*Period*Period*24.0*24.0*3600.0*3600.0/(4.0*np.pi*np.pi),1.0/3.0) )    

################################################################################

def Kepler(phi, ecen):
    phi=phi*180.0/np.pi
    for kk in range(len(phi)): 
        while(phi[kk]>360.0):
            phi[kk]=phi[kk]-360.0  
        while(phi[kk]<0.0):
            phi[kk]=phi[kk]+360.0       
        if(phi[kk]>180):  phi[kk]=float(phi[kk]-360)
        if(phi[kk]<-181.0 or phi[kk]>181.0):  
            print("Phi:  ",  phi[kk], ecen[kk])
            input("Enter a number ")
    phi=phi*np.pi/180.0##radian 
    ksi=phi; 
    for iw in range(Nb):
        term=2.0/(iw+1.0)*ss.jv( int(iw+1),(iw+1.0)*ecen)*np.sin((iw+1)*phi)
        if(iw==0):   term0=np.abs(term)
        ksi+=term
        if(np.mean(np.abs(term))<np.mean(abs(thre)*term0) and iw>5):  
            break 
    return(ksi)          
    
################################################################################

def ploto(x1, y1, z1, Vx, Delf1, Delf2, phase, Flux2, flag):       
    plt.cla()
    plt.clf()
    plt.figure(figsize=(8, 6))
    plt.plot(phase, x1, color="r",label="x1(LoS)",  lw=1.2, alpha=0.95)
    plt.plot(phase, y1, color="b",label="y1(RIGHT)",lw=1.2, alpha=0.95)
    plt.plot(phase, z1, color="g",label="z1(UP)",   lw=1.2, alpha=0.95)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.grid("True")
    plt.legend()
    fig=plt.gcf()
    fig.savefig("./exam/xyz{0:d}.jpg".format(flag), dpi=200)
    #############################################################
    plt.cla()
    plt.clf()
    plt.figure(figsize=(8,6))
    plt.plot(phase,Delf1, "r-",label=r"$Delf1$", lw=1.2, alpha=0.95)
    plt.plot(phase,Delf2, "b--",label=r"$Delf2$", lw=1.2, alpha=0.95)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.grid("True")
    plt.legend()
    fig=plt.gcf()
    fig.savefig("./exam/Delf{0:d}.jpg".format(flag), dpi=200)
    #############################################################
    plt.cla()
    plt.clf()
    plt.figure(figsize=(8, 6))
    plt.plot(phase, Vx/np.mean(np.abs(Vx)), color="r",label="Vx", lw=1.2, alpha=0.95)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.grid("True")
    fig=plt.gcf()
    fig.savefig("./exam/vx{0:d}.jpg".format(flag), dpi=200)
    #############################################################
    plt.cla()
    plt.clf()
    plt.figure(figsize=(8, 6))
    plt.errorbar(tim,Flux, yerr=errf, fmt=".",markersize=3.5,color='darkgreen',ecolor='g',elinewidth=0.5,capsize=0,alpha=0.5)
    plt.plot(Flux2[:,0],Flux2[:,1], "m-", label="median ",  lw=1.5, alpha=1.0)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.grid("True")
    fig=plt.gcf()
    fig.savefig("./exam/light{0:d}.jpg".format(flag), dpi=200)
    return(1)
 
################################################################################
def lnprob(p):  
    mbh, ecen, inc, teta, rstar, limb=p   
    if(mbh<0.0 or mbh>70.0 or ecen>0.85 or abs(inc)>5.0 or abs(teta)>180.0 or rstar<2.5 or rstar>4.0 or limb<0.35 or limb>0.65):
        return(-1.0)
    else: 
        return(+1.0)    
        
################################################################################
def sortarray(Astar0 , Flux2):
    l=np.argmax(Astar0)
    Tmod=tmod-tmod[l]
    for i in range(Nm): 
        if(Tmod[i]<-0.5): Tmod[i]=Tmod[i]+1.0
        if(Tmod[i]> 0.5): Tmod[i]=Tmod[i]-1.0
    ll=np.argsort(Tmod)
    for i in range(Nm): 
        l=int(ll[i])
        Flux2[i,0]=Tmod[l]
        Flux2[i,1]=Astar0[l]#fblend*Astar0[l]*(1.0+ Delf1[l]*Delf2[l])+1.0-fblend
    return(Flux2) 
################################################################################  
def orbit(PP, tim, Flux, errf, flag):
    MBH, ecen, inc, teta, Rstar, limb =PP  
    inc= float(inc*np.pi/180.0)
    teta=float(teta*np.pi/180.0)
    VBB.a1=limb;
    a= Thirdlaw(MBH)
    if(ecen<0.01): ksi=phi
    else:          ksi=Kepler(phi, ecen)
    x0=a*(np.cos(ksi)-ecen)
    y0=a*np.sin(ksi)*np.sqrt(1.0-ecen**2.0)
    y1=               y0*np.cos(teta)+x0*np.sin(teta)
    x1= np.cos(inc)*(-y0*np.sin(teta)+x0*np.cos(teta))
    z1=-np.sin(inc)*(-y0*np.sin(teta)+x0*np.cos(teta)) 
    ####################################
    ## Doppler Boosting 
    for s1 in range(Nm-1):   
        Vx[s1]= float(x1[s1+1]-x1[s1])/(tmod[s1+1]-tmod[s1])/(Period*24.0*3600.0)##[m/s]
    Vx[Nm-1]=Vx[Nm-2]+(Vx[Nm-2]-Vx[Nm-3])*(tmod[Nm-1]-tmod[Nm-2])/(tmod[Nm-2]-tmod[Nm-3])
    Delf1=Doppler(Vx)
    
    #####################################
    ## Ellipsoidal Variations
    dis= np.sqrt(x1**2.0 + y1**2.0 + z1**2.0)+1.0e-50
    phase=np.arccos(-x1/dis)      
    qm= MBH/Mstar
    cosi2= np.cos(inc)*np.cos(inc)
    li3= float(3.0-limb)
    ampl= Rstar*Rsun/dis
    gc= 0.01
    L0= abs(1.0+(15.0+limb)*(1.0+gc)*pow(ampl,3.0)*(2.0+5.0*qm)*(2.0-3.0*cosi2)/(60.0*(3.0-limb)) +9.0*(1.0-limb)*(3.0+gc)*pow(ampl,5.0)*qm*(8.0-40.0*cosi2+35.0*cosi2*cosi2)/(256.0*(3.0-limb)))+1.0e-50
    Delf2=((15*limb*(2.0+gc)*np.power(ampl,4.0)*qm*(4.0*np.cos(inc)-5.0*pow(np.cos(inc),3))/(32.0*li3))*np.cos(phase) + (-3.0*(15.0+limb)*(1.0+gc)*np.power(ampl,3.0)*qm*cosi2/(20.0*li3)-15.0*(1.0-limb)*(3.0+gc)*np.power(ampl,5.0)*qm*(6.0*cosi2-7.0*cosi2*cosi2)/(64.0*li3))*np.cos(2.0*phase)+ (-25.0*limb*(2.0+gc)*np.power(ampl,4.0)*qm*pow(np.cos(inc),3.0)/(32.0*li3))*np.cos(3.0*phase) + (105.0*(1.0-limb)*(3.0+gc)*np.power(ampl,5.0)*qm*pow(np.cos(inc),4.0)/(256.0*li3) )*np.cos(4.0*phase))/L0
      
    #####################################     
    # Self-Lensing  
    RE=const*np.sqrt(MBH)*np.sqrt(np.abs(x1)*Dl/(Dl+np.abs(x1))) + 1.0e-50#[m]
    rho= np.abs(Rstar*Rsun*Dl/(Dl+np.abs(x1))/RE)
    u=np.sqrt(y1**2.0 + z1**2.0)/RE
    for k in range(Nm): 
        if(x1[k]<0.0 and abs(phase[k])>np.pi/2.0 or RE[k]==0.0 or RE[k]<0.0 or a<=0.0 or rho[k]<0.0  or dis[k]<0.0 or dis[k]==0.0 or x1[k]>dis[k]): 
            print("Error x1, phase, kRE:  ",  x1[k], phase[k],  RE[k], limb)
            input("Enter a number  ")
        if(x1[k]<-0.001): 
            if(rho[k]>100.0):  
                if(u[k]<rho[k]): Astar[k]=float(1.0+2.0/rho[k]/rho[k])
                else:            Astar[k]=float(u[k]*u[k]+2.0)/np.sqrt(u[k]*u[k]*(u[k]*u[k]+4.0) )
            else:                Astar[k]=VBB.ESPLMag2(u[k],rho[k])    
        else:                    Astar[k]=1.0
              
              
    #####################################
    Flux2=np.zeros((Nm,2))
    Flux2=sortarray(Astar, Flux2) 
    chi2=0.0    
    for i in range(Nd):
        for j in range(Nm-1): 
            if(float((tim[i]-Flux2[j,0])*(tim[i]-Flux2[j+1,0]))<0.0 or tim[i]==Flux2[j,0]): 
                chi2+=float(Flux[i]-Flux2[j,1])**2.0/(errf[i]*errf[i])
                break
    if(flag>0): 
        out=ploto(x1/a, y1/a, z1/a, Vx*0.001, Delf1, Delf2, phase*180.0/np.pi, Flux2, flag)        
    return(chi2, Flux2)
################################################################################
def likelihood(p, tim, Flux, errf):
    test=lnprob(p)
    if(test<0.0): 
        return(-np.inf)
    else:
        chi2, Flux2=orbit(p,tim, Flux, errf,0)         
        return(-0.5*chi2)

################################################################################
def modelini(pp0):
    test= lnprob(pp0) 
    if(test>0.0):
        chi2,Flux2=orbit(pp0,tim, Flux, errf,0)          
    return(Flux2)    

################################  MAIN PROGRAM #################################
'''
##M_BH,    ecentricity,  inclination,  teta
## 0       1             2             3
# setting initial values  np.abs(np.random.normal(2.09774,0.05,Nwalkers))

#40.84175,  0.80966111,  1.36928,  0., 3.41675,   0.58572  
p0[:,0]=np.abs(np.random.normal(40.84,2.0,Nwalkers))## BH mass
p0[:,1]=np.abs(np.random.normal(0.81 ,0.2,Nwalkers))## ecentricity
p0[:,2]=       np.random.normal(1.37 ,0.5,Nwalkers)   ## inclination     
p0[:,3]=       np.random.normal(0.0  ,1.0,Nwalkers)   ## teta
p0[:,4]=np.abs(np.random.normal(3.42 ,0.1,Nwalkers))## Rstar
p0[:,5]=np.abs(np.random.normal(0.58 ,0.1,Nwalkers))## limb-darkening
#p0[:,0]=np.abs(np.random.normal(44.4,1.0,Nwalkers)) ## BH mass
#p0[:,1]=np.abs(np.random.normal(0.825,0.05,Nwalkers))    ## ecentricity
#p0[:,2]=np.random.normal(1.50,1.0,Nwalkers)      ## inclination     
#p0[:,3]=np.random.normal(-0.45,0.2,Nwalkers)   ## teta
#p0[:,4]=np.abs(np.random.normal(3.80,0.25,Nwalkers))## Rstar
#out=modelini(p0)
print("initla models are plotted")
sampler=emcee.EnsembleSampler(Nwalkers,NDIM,likelihood,args=(tim,Flux,errf),threads=4)
print("sampler is made !!!!!!")
xx=0
fil=open("./files/Chain1.dat", "a+")
fil.close()
for param, like, stat in sampler.sample(p0 ,iterations=Nstep, storechain=False):#, progress=True):
    xx+=1
    fil = open("./files/Chain1.dat","a+")
    ssa=np.concatenate((param.reshape((-1,NDIM)),like.reshape((-1,1))), axis=1)
    np.savetxt(fil,ssa,fmt ="%.5f   %.5f   %.5f   %.5f   %.5f   %.5f  %.5f")
    fil.close()
    if(xx%10==0): print("Step: " ,  xx,  np.mean(sampler.acceptance_fraction)   )
print("***************************** END OF MCMC **************************** ")
'''
################################################################################
col=["k", "b", "r", "g", "c", "y", "m"]
'''
fij=open("./exam/paramMBH.dat","w")
fij.close(); 
ns=22
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8, 6))
ax1=fig.add_subplot(111)
for i in range(ns): 
    mbh=float(0.7+(55.0-0.7)*i/ns/1.0)
    pp0=np.array([mbh,0.0,0.0,0.0,3.800000375436451,0.450000346354653])
    fij=open("./exam/paramMBH.dat","a+")
    np.savetxt(fij,pp0.reshape((-1,6)),fmt="%.5f  %.5f   %.5f  %.5f   %.5f   %.5f")
    fij.close(); 
    Flux2=modelini(pp0)
    if(i>6):  k=int(i%7)
    else:     k=i 
    if(i%2==0): plt.plot(Flux2[:,0]*Period*24.0, Flux2[:,1],color=col[k],ls="-", lw=float(0.14+i*0.08),alpha=1.0,label=r"$M_{\rm{c}}(M_{\odot})=$"+str(round(mbh,1)))
    else:       plt.plot(Flux2[:,0]*Period*24.0, Flux2[:,1],color=col[k],ls="-", lw=float(0.14+i*0.08),alpha=1.0)
    print("Mbh:  ", mbh)
plt.xlabel(r"$\rm{time}-t_{0}(\rm{hours})$", fontsize=18)
plt.ylabel(r"$\rm{Normalized}~\rm{Flux}$", fontsize=18, labelpad=0.01)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.xlim(-7.0, 7.0)
plt.title(r"$\epsilon=0,~i=0^{\circ},~\theta(\rm{deg})=0,~R_{\star}=3.8~R_{\odot},~\Gamma=0.45$",fontsize=18 ,color="k")
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=True)
plt.legend(prop={"size":16})
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./EffectMc.jpg",dpi=200)



################################################################################

fij=open("./exam/paramEcen.dat","w")
fij.close(); 
ns=22
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8, 6))
for i in range(ns): 
    ecen=float(0.+0.85*i/ns/1.0)
    pp0=np.array([40.00000364536,ecen,0.00003745235,0.0,3.800000375436451,0.450000346354653])
    fij=open("./exam/paramEcen.dat","a+")
    np.savetxt(fij,pp0.reshape((-1,6)),fmt="%.5f  %.5f   %.5f  %.5f   %.5f   %.5f")
    fij.close(); 
    Flux2=modelini(pp0)
    if(i>6):  k=int(i%7)
    else:     k=i 
    if(i%2==0):plt.plot(Flux2[:,0]*Period*24.0,Flux2[:,1],color=col[k],ls="-",lw=float(0.17+i*0.08),alpha=1.0,label=r"$\epsilon=$"+str(round(ecen,2)))
    else:      plt.plot(Flux2[:,0]*Period*24.0,Flux2[:,1],color=col[k],ls="-",lw=float(0.17+i*0.08),alpha=1.0)
    print("Ecen:  ", ecen)
plt.xlabel(r"$\rm{time}-t_{0}(\rm{hours})$", fontsize=18)
plt.ylabel(r"$\rm{Normalized}~\rm{Flux}$", fontsize=18, labelpad=0.01)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.xlim(-20.0 , 20.0)
plt.title(r"$M_{\rm{c}}=40~M_{\odot},~i=0^{\circ},~\theta(\rm{deg})=0,~R_{\star}=3.8~R_{\odot},~\Gamma=0.45$",fontsize=18,color="k")
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=True)
plt.legend(prop={"size":16})
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./EffectEcen.jpg",dpi=200)

################################################################################

fij=open("./exam/paramInc.dat","w")
fij.close(); 
ns=22
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8, 6))
for i in range(ns): 
    inc=float(0.0+3.0*i/ns/1.0)
    pp0=np.array([40.00000364536, 0.0, inc ,0.0, 3.800000375436451, 0.450000346354653])
    fij=open("./exam/paramInc.dat","a+")
    np.savetxt(fij,pp0.reshape((-1,6)),fmt="%.5f  %.5f   %.5f  %.5f   %.5f   %.5f")
    fij.close(); 
    Flux2=modelini(pp0)
    if(i>6):  k=int(i%7)
    else:     k=i 
    if(i%2==0):plt.plot(Flux2[:,0]*Period*24.0,Flux2[:,1],color=col[k],ls="-",lw=float(0.17+i*0.08),alpha=1.0,label=r"$i(\rm{deg})=$"+str(round(inc,2)))
    else:      plt.plot(Flux2[:,0]*Period*24.0,Flux2[:,1],color=col[k],ls="-",lw=float(0.17+i*0.08),alpha=1.0)
    print("inclination:  ", inc)
plt.xlabel(r"$\rm{time}-t_{0}(\rm{hours})$", fontsize=18)
plt.ylabel(r"$\rm{Normalized}~\rm{Flux}$", fontsize=18, labelpad=0.01)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.xlim(-5.0 , 5.0)
plt.title(r"$M_{\rm{c}}=40~M_{\odot},~\epsilon=0,~\theta(\rm{deg})=0,~R_{\star}=3.8~R_{\odot},~\Gamma=0.45$",fontsize=18,color="k")
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=True)
plt.legend(prop={"size":16})
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./EffectInc.jpg",dpi=200)



################################################################################

fij=open("./exam/paramRs.dat","w")
fij.close(); 
ns=22
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8, 6))
for i in range(ns): 
    Rs=float(3.4+0.5*i/ns/1.0)
    pp0=np.array([40.00000364536, 0.0, 0.0 ,0.0, Rs, 0.450000346354653])
    fij=open("./exam/paramRs.dat","a+")
    np.savetxt(fij,pp0.reshape((-1,6)),fmt="%.5f  %.5f   %.5f  %.5f   %.5f   %.5f")
    fij.close(); 
    Flux2=modelini(pp0)
    if(i>6):  k=int(i%7)
    else:     k=i 
    if(i%2==0):plt.plot(Flux2[:,0]*Period*24.0,Flux2[:,1],color=col[k],ls="-",lw=float(0.17+i*0.08),alpha=1.0,label=r"$R_{\star}(R_{\odot})=$"+str(round(Rs,2)))
    else:      plt.plot(Flux2[:,0]*Period*24.0,Flux2[:,1],color=col[k],ls="-",lw=float(0.17+i*0.08),alpha=1.0)
    print("Rs:  ", Rs)
plt.xlabel(r"$\rm{time}-t_{0}(\rm{hours})$", fontsize=18)
plt.ylabel(r"$\rm{Normalized}~\rm{Flux}$", fontsize=18, labelpad=0.01)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.xlim(-5.0,  5.0)
plt.title(r"$M_{\rm{c}}=40~M_{\odot},~\epsilon=0,~i=0^{\circ},~\theta(\rm{deg})=0,~\Gamma=0.45$",fontsize=18,color="k")
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=True)
plt.legend(prop={"size":16})
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./EffectRS.jpg",dpi=200)
'''

################################################################################
'''
fij=open("./exam/paramTET.dat","w")
fij.close(); 
ns=22
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8, 6))
for i in range(ns): 
    tet=float(179.63452462523746*i/ns/1.0)
    pp0=np.array([40.00364536, 0.6549563475645,0.00003645364,tet, 3.8000375436451, 0.44475647563])
    fij=open("./exam/paramTET.dat","a+")
    np.savetxt(fij,pp0.reshape((-1,6)),fmt="%.5f  %.5f   %.5f  %.5f   %.5f   %.5f")
    fij.close(); 
    Flux2=modelini(pp0)
    
    for k in range(Nm): 
        if(abs(Flux2[k,1]-Flux2[k-1,1])>abs(Flux2[k-1,1]-Flux2[k-2,1]) and abs(Flux2[k,0])<float(0.9/Nm/1.0) ):  
            print(Flux2[k,1], Flux2[k-1,1], Flux2[k-2,1], Flux2[k,0])
            Flux2[k,1]=float(Flux2[k-1,1]+Flux2[k-2,1]+Flux2[k-3,1])/3.0

    
    if(i>6):  k=int(i%7)
    else:     k=i 
    if(i%2==0):plt.plot(Flux2[:,0]*Period*24.0,Flux2[:,1],color=col[k],ls="-",lw=float(0.17+i*0.08),alpha=1.0,label=r"$\theta(\rm{deg})=$"+str(int(tet)))
    else:      plt.plot(Flux2[:,0]*Period*24.0,Flux2[:,1],color=col[k],ls="-",lw=float(0.17+i*0.08),alpha=1.0)
    print("TETA:  ", tet)
plt.xlabel(r"$\rm{time}-t_{0}(\rm{hours})$", fontsize=18)
plt.ylabel(r"$\rm{Normalized}~\rm{Flux}$", fontsize=18, labelpad=0.01)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.xlim(-10.0,  10.0)
plt.title(r"$M_{\rm{c}}=40~M_{\odot},~\epsilon=0.65,~i=0^{\circ},~R_{\star}(R_{\odot})=3.8,~\Gamma=0.45$",fontsize=18,color="k")
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=True)
plt.legend(prop={"size":14})
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./EffectTET.jpg",dpi=200)

'''

################################################################################

fij=open("./exam/paramlimb.dat","w")
fij.close(); 
ns=22
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8, 6))
for i in range(ns): 
    limb=float(0.38+0.17*i/ns)
    pp0=np.array([40.00364536, 0.0, 0.00, 0.0 , 3.8000375436451, limb])
    fij=open("./exam/paramlimb.dat","a+")
    np.savetxt(fij,pp0.reshape((-1,6)),fmt="%.5f  %.5f   %.5f  %.5f   %.5f   %.5f")
    fij.close(); 
    Flux2=modelini(pp0)
    
    for k in range(Nm): 
        if(abs(Flux2[k,1]-Flux2[k-1,1])>abs(Flux2[k-1,1]-Flux2[k-2,1]) and abs(Flux2[k,0])<float(0.9/Nm/1.0) ):  
            print(Flux2[k,1], Flux2[k-1,1], Flux2[k-2,1], Flux2[k,0])
            Flux2[k,1]=float(Flux2[k-1,1]+Flux2[k-2,1]+Flux2[k-3,1])/3.0
           
    if(i>6):  k=int(i%7)
    else:     k=i 
    if(i%2==0):plt.plot(Flux2[:,0]*Period*24.0,Flux2[:,1],color=col[k],ls="-",lw=float(0.17+i*0.08),alpha=1.0,label=r"$\Gamma=$"+str(round(limb,2)))
    else:      plt.plot(Flux2[:,0]*Period*24.0,Flux2[:,1],color=col[k],ls="-",lw=float(0.17+i*0.08),alpha=1.0)
    print("limb:  ", limb)
plt.xlabel(r"$\rm{time}-t_{0}(\rm{hours})$", fontsize=18)
plt.ylabel(r"$\rm{Normalized}~\rm{Flux}$", fontsize=18, labelpad=0.01)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.xlim(-5.0,  5.0)
plt.title(r"$M_{\rm{c}}=40~M_{\odot},~\epsilon=0,~i=0^{\circ},~\theta=0^{\circ},~R_{\star}(R_{\odot})=3.8$",fontsize=18,color="k")
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=True)
plt.legend(prop={"size":14})
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./Effectlimb.jpg",dpi=200)



################################################################################




