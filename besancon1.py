import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import rcParams
import warnings
warnings.filterwarnings("ignore")
rcParams["font.size"] = 13
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
##      0     1      2        3       4      5     6   7     8     9      10   11     12      13      14       15        16
##     V     B-V     U-B     V-I     V-K     Mv   CL Typ    Teff   logg  Pop  Age     Mass   Mbol    Radius   [M/H]    [a/Fe]  longitude   latitude   RAJ2000     DECJ2000   Dist     x_Gal    y_Gal    z_Gal     Av       errMv     errMass     errAge     errTeff    errLogg    errMet     errAlphaFe     errBand_V  errBand_B  errBand_U  errBand_I  errBand_K  

#     V     B-V     U-B     V-I     V-K     Mv   CL Typ    Teff   logg  Pop  Age     Mass   Mbol    Radius   [M/H]    [a/Fe]  longitude   latitude   RAJ2000     DECJ2000   Dist     x_Gal    y_Gal    z_Gal     Av       errMv     errMass     errAge     errTeff    errLogg    errMet     errAlphaFe     errBand_V  errBand_B  errBand_U  errBand_I  errBand_K 

#19.488  1.541   1.165   2.424   4.709  11.99  5 7.30    3321.  4.86   2  0.1072  0.543    9.798    0.287   0.086   0.017     0.406932   -4.695548  271.325409  -30.946587   0.3123  -7.6888   0.0022  -0.0106   0.043    0.0000     0.1000     0.1000    50.0000     0.1000     0.1000     0.0200    -0.0183    -0.0105    -0.0042    -0.0006     0.0183

'''
fij=open("./Besanconx.dat","a+")
fij.close(); 
f0=open("./outBosancon5.txt","r")
nr=sum(1 for line in f0)
par=np.zeros(( nr , 38))
par=np.loadtxt("./outBosancon5.txt") 
G=   6.67430*pow(10.0,-11.0)
Msun=1.989*pow(10.0,30.0)
Rsun =6.9634*pow(10.0,8.0)
nn=np.zeros((9))
good=0
for i in range(nr):  
    mass, age, metal1, metal2, logg, teff,Mv,radius,cl=par[i,12], par[i,11], par[i,15], par[i,16],  par[i,9], par[i,8], par[i,5], par[i,14], par[i,6]
    if(mass >1.5 and mass<2.5 and radius<5.0):  
        good+=1
        nn=np.array([mass, age, metal1, metal2, logg, teff,Mv,radius, cl])
        #GG=pow(10.0,logg)*0.01
        #gra= G*mass*Msun*(radius*Rsun)**-2.0
        #print ("Gravity:  ",  GG,    gra  )
        fij=open("./Besanconx.dat","a+")
        np.savetxt(fij, nn.reshape((1,9)),fmt= "%.3f   %.4f    %.3f     %.3f   %.3f    %.1f    %.2f    %.4f  %d")
        fij.close()
print (good)        
input("Enter a number ")
'''

f0=open("./Besanconx.dat","r")
num0=sum(1 for line in f0)
par=np.zeros(( num0, 9))
par=np.loadtxt("./Besanconx.dat") 
file3=open("Besancon_sort.dat","w")
file3.close()
massi=np.zeros((num0)); 
csort=np.zeros((num0))
massi=par[:,0] 
csort=np.argsort(massi); 
for i in range(num0):
    if(par[int(csort[i]),7]>2.0): ### limitation of stellar radius 
        file3=open("Besancon_sort.dat","a") 
        np.savetxt(file3,par[int(csort[i]),:].reshape(1,9),fmt='%.3f   %.4f    %.3f     %.3f   %.3f    %.1f    %.2f    %.4f  %d') 
        file3.close()


plt.cla()
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(par[:,0], par[:,1], "ro",label=r"$Age-Magg$", alpha=0.95)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
#plt.yscale('log')
plt.grid("True")
plt.xlabel(r"$\rm{Mass}(M_{\odot})$", fontsize=18)
plt.xlabel(r"$\rm{Age}(\rm{Gyrs})$", fontsize=18)
plt.grid(linestyle='dashed')
plt.legend()
fig3=plt.gcf()
fig3.savefig("./MassAgeB.jpg", dpi=200)

input("Enter a number ")


f1=open("Besancon_sort.dat","r")
num0=sum(1 for line in f1)
################################################################################
par=np.zeros(( num0, 9))
par=np.loadtxt("./Besancon_sort.dat")     
file4=open("Besancon_sort2.dat","w")
file4.close()
massr=np.arange(1.5, 2.5, 0.001, dtype=float)
count= np.zeros((len(massr)))   
mass0= np.zeros((len(massr))) 
kk=0
print (massr,  len(massr))
arry=np.zeros((3000,9))
for i in range(len(massr)):
    k=0
    for j in range(num0): 
        if(abs(par[j,0]-massr[i])<0.0008): 
            arry[k,:]=par[j,:]
            k+=1
            kk+=1
            mass0[i]=par[j,0]
    #print ("mass:  ",  massr[i],   k)        
    
    #if(i==0):  count[i]=0
    #else:      
    count[i]=kk # + count[i-1]

    
    
    csort=np.zeros((k))        
    csort=np.argsort(arry[:k,1]) 
    for l in range(k):  
        file4=open("Besancon_sort2.dat","a+")
        np.savetxt(file4,arry[int(csort[l]),:].reshape(1,9),fmt='%.3f   %.4f    %.3f     %.3f   %.3f    %.1f    %.2f    %.4f %d') 
        file4.close()
################################################################################    

print("Count:  ",  count)
print (len(count))    
file4=open("massC.dat","w")
for i in range(len(count)): 
    if(mass0[i]>0.0): 
        tt=np.array([ mass0[i], count[i] ])
        np.savetxt(file4,tt.reshape(1, 2),fmt='%.3f    %d')
#np.savetxt(file4,mass0.reshape(1, len(massr)),fmt='%.3f, ') 
file4.close()    
input("Enter a number ")

##mass, age, metal1, metal2, logg, teff,Mv,radius
f0=open("./Besancon_sort2.dat","r")
num0=sum(1 for line in f0)
par=np.zeros(( num0, 9))
par=np.loadtxt("./Besancon_sort2.dat") 

nl= 9586
limb=np.zeros(( nl , 7))
limb=np.loadtxt("./Limb_Linear_Tess24.dat") 
file4=open("Besancon_sort2_Limb.dat","w")
file4.close()
met=[-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.5, 1.0]
## logg  Teff   Z   LHP   uLSM    uPCM   chi2   Type   Mod 
##  0     1     2    3    4       5      6       7     8 



for i in range(num0): 
    j0=0
    for j in range(len(met)-1):
        if( float((par[i,2]-met[j])*(par[i,2]-met[j+1]))<0.0  or par[i,2]==met[j] ): 
            j0=int(j)
            break
    #print ("metal:  ", par[i,2],  met[j0] )
    #input("Enter a number ")        
    array=np.zeros((2500, 7))
    k=0;         
    for j in range(nl): 
        if(abs(met[j0]-limb[j,2])<0.05): 
            array[k, :]= limb[j,:]       
            k+=1
    #print ("Number of elements:  ",   k)        
    dism=10000054;  
    k0=0; 
    for j in range(k):
        dis=np.sqrt( (par[i,5]- array[j,1])**2.0 + (par[i,4]-array[j,0])**2.0 ) 
        if(dis<dism): 
            k0=j
            dism=dis                      
    #print ("Besancon: Teff, logg, metal  ", par[i,5],   par[i,4],    par[i,2] )
    #print ("Limb_TESS:Teff, logg, metal  ", array[k0,1],  array[k0,0], array[k0,2])
    #input("Enter a number ")
    mass, radius, teff=  par[i,0], par[i,7], par[i,5]
    if(abs(mass-2.08)<0.06 and radius>3.39 and radius<3.89 and teff>7455.00 and teff<7986.99):  
        print ( par[i,0],  par[i,1], par[i,2], par[i,4],  par[i,5], par[i,6], par[i,7], array[k0,4] , par[i,8] )
    
    if(par[i,7]>1.0 and par[i,7]<5.0 and par[i,0]>1.5 and par[i,0]<2.5):
        file4=open("Besancon_sort2_Limb.dat","a+")
        save= np.array([ par[i,0],  par[i,1], par[i,2], par[i,4],  par[i,5], par[i,6], par[i,7], array[k0,4] , par[i,8]])
        np.savetxt(file4,save.reshape((1, 9)), fmt='%.3f   %.4f    %.3f    %.3f    %.1f    %.2f    %.4f   %.4f  %d')
        file4.close() 



##mass, age, metal1, metal2, logg, teff,Mv,radius
























