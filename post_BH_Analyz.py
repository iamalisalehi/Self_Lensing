from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import rcParams 
import matplotlib as mpl
import scipy.special as ss
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import emcee, corner
from matplotlib import gridspec
import warnings
from matplotlib import cm
warnings.filterwarnings("ignore")
rcParams["font.size"] = 13
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
from numpy import matrix
import VBBinaryLensingLibrary as vb
VBB=vb.VBBinaryLensing()
VBB.Tol=1.0e-3; 
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
def ploto(x1, y1, z1, Vx, Delf1, Delf2, phase, flag):       
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
    fig.savefig("./Figs/xyz{0:d}.jpg".format(flag), dpi=200)
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
    fig.savefig("./Figs/Delf{0:d}.jpg".format(flag), dpi=200)
    #############################################################
    plt.cla()
    plt.clf()
    plt.figure(figsize=(8, 6))
    plt.plot(phase, Vx/np.mean(np.abs(Vx)), color="r",label="Vx", lw=1.2, alpha=0.95)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.grid("True")
    fig=plt.gcf()
    fig.savefig("./Figs/vx{0:d}.jpg".format(flag), dpi=200)
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
        out=ploto(x1/a, y1/a, z1/a, Vx*0.001, Delf1, Delf2, phase*180.0/np.pi, flag)        
    return(chi2, Flux2)
################################################################################
def likelihood(p, tim, Flux, errf):
    test=lnprob(p)
    if(test<0.0): 
        return(-np.inf)
    else:
        chi2, Flux2=orbit(p,tim, Flux, errf,0)         
        return(-1.0*chi2)
################################  MAIN PROGRAM #################################

a= np.zeros((NDIM,3))#median
c= np.zeros(NDIM)## MBH,   Ecentricity, inc, teta, Rstar, Gamma
fdd=open("./files/chaintot.dat","r")
nns=sum(1 for line in fdd)
Chain= np.zeros((nns,NDIM+1))
Chain=np.loadtxt("./files/chaintot.dat")  
samples1=Chain[:,:NDIM].reshape((-1,NDIM))
Fluxm= np.zeros((Nm,2));
a=map(lambda v:(v[1],v[2]-v[1],v[1]-v[0]),zip(*np.percentile(samples1, [16, 50, 84],axis=0)))
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
print("chi2 of median  model: ",    abs(chi2c),  asemic)
print("circular velocity(km/s):  ",2.0*np.pi*asemic*AU*0.001/Period/(24.0*3600.0))
print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")

################################################################################
b=np.zeros(NDIM)
tt=np.argmax(Chain[:,NDIM])
b=Chain[tt,:NDIM]
b=np.array([43.38182, 0.82470, 1.53584, -0.52086, 3.8, 0.45])
Fluxb=np.zeros((Nm,2))
chi2b,Fluxb=orbit(b, tim, Flux, errf ,2)
asemib=Thirdlaw(b[0])/AU
save2=open("./files/BestsN.txt","a+")
np.savetxt(save2,b,delimiter='  ',fmt='%.10f',comments='#',header='The best fitted (Chi2_min) parameters')
save2.close()
print("*the best model : " ,  b) 
print("Chi2 of the best model:  "     , abs(chi2b),   asemib)
print("Circular velocity(km/s):",2.0*np.pi*asemib*AU*0.001/Period/(24.0*3600.0))
print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")

################################################################################
    
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8, 6))
spec2=gridspec.GridSpec(ncols=1, nrows=1)
ax =fig.add_subplot(spec2[0,0])
ax.errorbar(tim,Flux, yerr=errf, fmt=".",markersize=4.9,color='darkgreen',ecolor='g',elinewidth=0.6,capsize=0,alpha=0.5)
ax.plot(Fluxb[:,0],Fluxb[:,1],"m-",label=r"$\rm{MCMC}_{A}$", lw=1.3, alpha=1.0)
ax.plot(Fluxm[:,0],Fluxm[:,1],"k-",label=r"$\rm{MCMC}_{B}$", lw=1.3, alpha=1.0)
plt.xlabel(r"$\rm{time}/Period$", fontsize=18)
plt.ylabel(r"$\rm{Normalized}~\rm{Flux}$", fontsize=18, labelpad=0.01)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.xlim(-0.5,0.5)
plt.ylim(0.9985,1.0063)
plt.legend()
plt.legend(prop={"size":17.0})
left, bottom, width, height = [0.12, 0.53, 0.4, 0.45]
axins = fig.add_axes([left, bottom, width, height])
plt.setp(axins.get_xticklabels(), visible=False)
plt.setp(axins.get_yticklabels(), visible=False)
axins.yaxis.get_major_locator().set_params(nbins=7)
axins.xaxis.get_major_locator().set_params(nbins=7)
axins.errorbar(tim,Flux, yerr=errf, fmt=".",markersize=4.9,color='darkgreen',ecolor='g',elinewidth=0.6,capsize=0,alpha=0.5)
axins.plot(Fluxb[:,0],Fluxb[:,1], "m-",  lw=1.3, alpha=1.0)
axins.plot(Fluxm[:,0],Fluxm[:,1], "k-",  lw=1.3, alpha=1.0)
axins.set_ylim(0.998,1.0065)
axins.set_xlim(-0.1,0.1)
axins.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
axins.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
fig3=plt.gcf()
fig3.tight_layout(pad=0.2)
fig3.savefig("./Figs/BestModeln.jpg", dpi=200)
print(">>>>>>>>>>>>>>>> The light curve was made <<<<<<<<<<<<<<<<<<<<<<<<")

################################################################################

nam= [r"$M_{\rm{BH}}(M_{\odot})$",r"$\epsilon$",r"$|i^{\circ}|$",r"$|\theta^{\circ}|$", r"$R_{\star}(R_{\odot})$", r"$\Gamma$"]
samples1[:,2]=np.abs(samples1[:,2])
samples1[:,3]=np.abs(samples1[:,3])
plt.cla()
plt.clf()
plt.figure(figsize=(6,6))
'''
fig=corner.corner(samples1,labels=nam,quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 23},color="k",truths=b, truth_color="b")
fig.gca().annotate("",xy=(1.0, 1.0), xycoords="figure fraction", xytext=(-20,-10), textcoords="offset points", ha="right", va="top")
fig.subplots_adjust(right=1.5,top=1.5)
for ax in fig.get_axes():
    ax.tick_params(axis='both', labelsize=25)
plt.gcf()
fig.savefig("./Figs/corner2n.jpg",dpi=200, pad_inches=0.3,bbox_inches='tight')
print(">>>>>>>>corner is plotted <<<<<<<<<<<<<<<<<<")
'''
fig = corner.corner(samples1[:,:int(NDIM-1)], show_titles=True,title_kwargs={"fontsize":22},range=[(37.0, 45.0), (0.78, 0.84), (1.3, 1.6), (0, 0.9),(3.4,3.85)], labels=nam[:int(NDIM-1)],label_kwargs={"fontsize":22}, ticks_kwargs={"fontsize":24},quantiles=[0.16,0.5,0.84],plot_datapoints=True)##, truths=b,truth_color="m")
fig.savefig("./Figs/trianglen.jpg", dpi=200)
print(">>>>>>>>>>>>>>>> The triangle plot was made <<<<<<<<<<<<<<<<<<<<<<<<")
################################################################################

unique, index=np.unique(-2.0*Chain[:,NDIM],axis=0, return_index=True)
arra=np.zeros((len(index),NDIM+1))
print("total: ",  len(Chain[:,NDIM]),  "New not-repeated ones: ",   len(index) )
for i in range(len(index)):  
    arra[i,:]=Chain[int(index[i]), :]
    arra[i,NDIM]=abs(arra[i,NDIM])
    
lk=np.argsort(arra[:,NDIM])    
print(lk)
for i in range(10): 
    print("Sorted for chi2:  ", arra[int(lk[i]),:], "Chi2: ",    round(arra[lk[i],NDIM],1),    int(lk[i])) 
    print("************************************************")   
#input("Enter a number ")    
    
test=open("./files/tops.txt","w")
test.close()
nl=int(250)
Fluxt=np.zeros((Nm,2, nl))
chi2m=np.zeros((nl))
for i in range(nl): 
    bt=arra[int(lk[i*10]),:NDIM]
    chi2t,Fluxt[:,:,i]=orbit(bt, tim, Flux, errf ,0)
    chi2m[i]=abs(chi2t)
    asemit=Thirdlaw(bt[0])/AU
    print("parameters:  ",    abs(arra[int(lk[i*10]),NDIM]*2.0),    float(chi2t),  asemit)
    print("***************************************************************")   
    test=open("./files/tops.txt","a+")
    sav=np.array([i*4, bt[0],  bt[1],   bt[2],  bt[3],  bt[4],  abs(arra[int(lk[i*10]),NDIM]*2.0)   ])
    np.savetxt(test, sav.reshape((-1,7)), fmt= "%d   %.5f    %.5f    %.5f   %.5f   %.5f   %.2f")
    test.close()
    
mic=np.min(chi2m)
leg=np.max(chi2m)-mic
ticp=np.array([0.1,  0.3, 0.5,  0.7, 0.9  ])
ticl=np.round(ticp*leg+mic,1)
#for i in range(5): 
#    ticl[i]=int(ticl[i])
print ("tics_position:  ", ticp)
print ("tics_label:  :  ", ticl)
cmap = mpl.cm.get_cmap('seismic')
rgba =cmap(np.linspace(0, 1.0, len(chi2m)))
plt.cla()
plt.clf()
fig, ax = plt.subplots(figsize=(8, 6))
for j in range(nl): 
    ax.plot(Fluxt[:,0,j],Fluxt[:,1,j], color=cmap((chi2m[j]-mic)/leg), lw=0.9, alpha=1.0)
plt.errorbar(tim,Flux,yerr=errf, fmt=".",markersize=8.485,color='g',ecolor='#C1FFC1',elinewidth=1.85,capsize=0,alpha=0.45)    
ssww=mpl.cm.ScalarMappable(cmap=mpl.colors.ListedColormap(rgba))
ssww.set_array([])
cbar=plt.colorbar(ssww,orientation='vertical',ticks=ticp, pad=0.01)
cbar.set_ticklabels(ticl)
cbar.ax.tick_params(labelsize=16) 
cbar.set_label(label=r"$\chi^{2}$",size=18,weight='bold')
plt.xlabel(r"$\rm{time}/Period$", fontsize=18)
plt.ylabel(r"$\rm{Normalized}~\rm{Flux}$", fontsize=18, labelpad=0.01)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.ylim(0.998,1.0065)
plt.xlim(-0.13, 0.13)
plt.legend()
plt.legend(prop={"size":17.0})
fig3=plt.gcf()
fig3.tight_layout(pad=0.2)
fig3.savefig("./Figs/top200n.jpg", dpi=200)
































################################################################################



























