from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import rcParams 
import matplotlib as mpl
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('seaborn-poster')
import scipy.special as ss
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import emcee, corner
from matplotlib import gridspec
from scipy.ndimage.filters import gaussian_filter
import warnings
from matplotlib import cm
warnings.filterwarnings("ignore")
rcParams["font.size"] = 13
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
mpl.rcParams['legend.fontsize'] = 10

import emcee, corner
from numpy import matrix
import VBBinaryLensingLibrary as vb
VBB=vb.VBBinaryLensing()
VBB.Tol=1.0e-3; ### accuracy
VBB.SetLDprofile(VBB.LDlinear);
VBB.LoadESPLTable("./ESPL.tbl"); 
################################################################################
###        Constant            ###
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




Mstar=2.08
Tstar=7542.0
fblend=1.0
Period=11.393423143561519#[days]
Dl=float(1.0/0.978796586561136)*KPC## meter
tp=654.0*Period/3000.0
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
tim=tim/Period
thre=float(0.001)
Nb=1000
################################################################################
NDIM=4+2
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
def ploto(x1, y1, z1, timm, fluxx):
    for i in range(50):
        j=int(Nm/50*i)        
        plt.cla()
        plt.clf()
        fig = plt.figure(figsize=plt.figaspect(0.5))
        ax1 =fig.add_subplot(1,2,1, projection='3d')
        ax1.grid()
        ax1.plot(x1, y1, z1, color="g",label=r"$\rm{Stellar}~\rm{Orbit}$",lw=1.6, alpha=1.0)
        ax1.scatter(0.0, 0.0, 0.0, color="k", s=50, label=r"$\rm{Compact}~\rm{Object}$")
        ax1.scatter(x1[j],y1[j],z1[j],marker= "*",facecolors='r', edgecolors='k',s=100,label=r"$\rm{Source}~\rm{star}$")
        ax1.set_xlabel(r"$x_{\rm{o}}/a$", labelpad=22, fontsize=18)
        ax1.set_ylabel(r"$y_{\rm{o}}/a$", labelpad=17, fontsize=18)
        ax1.set_zlabel(r"$z_{\rm{o}}/a$", labelpad=12, fontsize=18)
        ax1.legend(prop={"size":14.0})
        ax1.view_init(azim=45+15, elev=30)
        #fig=plt.gcf()
        #fig.tight_layout()
        #fig.savefig("./Figs/orbit{0:d}.jpg".format(i), dpi=200)
        #############################################################
        #plt.cla()
        #plt.clf()
        #plt.figure(figsize=(8,6))
        #fig=plt.figure(figsize=(8, 6))
        #fig,ax=plt.subplots(figsize=(8,6))
        ax2=fig.add_subplot(1,2,2)
        ax2.plot(timm*Period,fluxx,color="k", lw=1.2, alpha=1.0)
        ax2.errorbar(tim*Period,Flux,yerr=errf,fmt=".",markersize=5.4,color='g',ecolor='darkgreen',elinewidth=0.5, capsize=0,alpha=0.35) 
        ax2.scatter(timm[j]*Period,fluxx[j],marker= "o",facecolors='r', edgecolors='k',s= 60.0,alpha=1.0)
        ax2.set_xlabel(r"$\rm{time}-t_{0}(\rm{days})$", fontsize=18)
        ax2.set_ylabel(r"$\rm{Normalized}~\rm{Flux}$", fontsize=18, labelpad=0.01)
        plt.xticks(fontsize=17, rotation=0)
        plt.yticks(fontsize=17, rotation=0)
        ax2.set_ylim(0.998,1.0065)
        ax2.set_xlim(-Period/2.0, Period/2.0)
        ax2.legend()
        ax2.legend(prop={"size":17.0})
        fig=plt.gcf()
        fig.tight_layout()
        fig.savefig("./Figs/MO{0:d}.jpg".format(i),dpi=200)
        print(i, "******************************************")
    return(1)
################################################################################
def lnprob(p):  
    mbh, ecen, inc, teta, rstar=p   
    if(mbh<2.0 or mbh>60.0 or ecen>0.85 or inc>10.0 or teta>180.0 or  teta<0.0 or inc<0.0 or rstar<1.5 or rstar>4.5):
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
        Flux2[i,1]=Astar0[l]
    return(Flux2) 
################################################################################  
def orbit(PP, tim, Flux, errf, flag):
    MBH, ecen, inc, teta, Rstar, limb =PP  
    inc= float(inc*np.pi/180.0)
    teta=float(teta*np.pi/180.0)
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
    VBB.a1=limb;
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
        out=ploto(x1/a, y1/a, z1/a, tmod, Astar )        
    return(chi2, Flux2)
################################################################################
def likelihood(p, tim, Flux, errf):
    test=lnprob(p)
    if(test<0.0): 
        return(-np.inf)
    else:
        chi2, Flux2=orbit(p,tim, Flux, errf,0)         
        return(-0.5*chi2)
################################ MAIN PROGRAM #################################
a= np.zeros((NDIM,3))#median
c= np.zeros(NDIM)## MBH,   Ecentricity, inc, teta, Rstar
fdd=open("./files/chaintot.dat","r")
nns=sum(1 for line in fdd)
Chain= np.zeros((nns,NDIM+1))
Chain=np.loadtxt("./files/chaintot.dat") 
samples1=Chain[:,:NDIM].reshape((-1,NDIM))
Fluxm= np.zeros((Nm,2));
a=map(lambda v:(v[1],v[2]-v[1],v[1]-v[0]),zip(*np.percentile(samples1, [16, 50, 84],axis=0)))
#print(zip(*np.percentile(samples1, [16, 50, 84],axis=0)) )
c=matrix(a).transpose()[0].getA()[0]


c[0]=c[0]-1.5-0.1
c[4]=c[4]-0.01
c[1]=c[1]-0.01


chi2c, Fluxm=orbit(c, tim, Flux, errf ,1)
asemic= Thirdlaw(c[0])/AU
save2=open("./files/BestsN.txt","w")
np.savetxt(save2,a,delimiter=' ',fmt='%.10f',comments='#',header='The best fitted (median) parameters')
save2.close()
for i in range(NDIM): print(c[i], "   +  ", a[i][1], "  -  ", a[i][2] )
print("parameters of the median model :  ",     c )
print("chi2 of median  model: ",    chi2c/2.0,  asemic)
print("circular velocity(km/s):  ",2.0*np.pi*asemic*AU*0.001/Period/(24.0*3600.0))
print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
#b=np.zeros(NDIM)
#tt=np.argmax(Chain[:,NDIM])
#b=Chain[tt,:NDIM]
#Fluxb=np.zeros((Nm,2))
#chi2b,Fluxb=orbit(b, tim, Flux, errf ,0)
#asemib=Thirdlaw(b[0])/AU
#save2=open("./files/BestsN.txt","a+")
#np.savetxt(save2,b,delimiter='  ',fmt='%.10f',comments='#',header='The best fitted (Chi2_min) parameters')
#save2.close()
#print("*the best model : " ,  b) 
#print("Chi2 of the best model:  "     , chi2b/2.0,   asemib)
#print("Circular velocity(km/s):",2.0*np.pi*asemib*AU*0.001/Period/(24.0*3600.0))
print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
################################################################################    

























