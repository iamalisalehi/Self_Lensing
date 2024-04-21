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
VBB.Tol=1.0e-3; ### accuracy
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

NDIM=4+1
Nwalkers=24
Nstep= 300000
p0=np.zeros((Nwalkers,NDIM))
Nm=int(3000)
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
    tmod[i]=float(-0.5+i/Nm/1.0) #[-0.5 , 0.5]
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
    mbh, ecen, inc, teta, rstar=p 
    limb=0.45  
    if(mbh<20.0 or mbh>50.0 or ecen>0.84 or abs(inc)>3.5 or abs(teta)>3.0 or rstar<3.38 or rstar>3.82 or limb<0.4 or limb>0.55):
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
    MBH, ecen, inc, teta, Rstar=PP  
    inc= float(inc*np.pi/180.0)
    teta=float(teta*np.pi/180.0)
    limb=0.45
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
def modelini(p0):
    fij=open("./exam/param2.dat","w")
    fij.close(); 
    for it in range(Nwalkers):
        print("initial parameters:  ",  p0[it,0], p0[it,1],  p0[it,2], p0[it,3] ) 
        test= lnprob(p0[it,:]) 
        if(test>0.0):
            chi2, Flux2=orbit(p0[it,:],tim, Flux, errf,it+Nwalkers)          
    return(1)    

################################  MAIN PROGRAM #################################
##M_BH,    ecentricity,  inclination,  teta, Rstar
## 0       1             2             3      4
#40.84175,  0.80966111,  1.36928,  0., 3.41675,   0.58572 
 
p0[:,0]=np.abs(np.random.rand(Nwalkers)*15.0+33.0)## BH mass     [33.0,48.0]
p0[:,1]=np.abs(np.random.rand(Nwalkers)*0.3+0.55) ## ecentricity [0.55,0.85]
p0[:,2]=np.abs(np.random.rand(Nwalkers)*1.7+0.50) ## inclination [0.5,2.2]    
p0[:,3]=np.abs(np.random.rand(Nwalkers)*2.5+0.00) ## teta        [0.0,2.5]
p0[:,4]=np.abs(np.random.rand(Nwalkers)*0.4+3.40) ## Rstar       [3.4,3.8]

#out=modelini(p0)
print("initla models are plotted")
sampler=emcee.EnsembleSampler(Nwalkers,NDIM,likelihood,args=(tim,Flux,errf),threads=24)
print("sampler is made !!!!!!")
xx=0
fil=open("./files/Chain1.dat", "w")
fil.close()
for param, like, stat in sampler.sample(p0 ,iterations=Nstep, storechain=False):#, progress=True):
    xx+=1
    fil = open("./files/Chain1.dat","a+")    
    ssa=np.concatenate(( param.reshape((-1,NDIM)) , like.reshape((-1,1)) ),axis=1)
    np.savetxt(fil,ssa,fmt ="%.5f   %.5f   %.5f   %.5f   %.5f  %.5f")
    fil.close()
    if(xx%10==0): print("Step: " ,  xx,  np.mean(sampler.acceptance_fraction))
print("***************************** END OF MCMC **************************** ")

################################################################################













