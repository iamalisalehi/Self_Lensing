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

fij=open("./Besancon75.dat","a+")
fij.close(); 



f0=open("./out_last.dat","r")
nr=sum(1 for line in f0)
par=np.zeros(( nr , 38))
par=np.loadtxt("./out_last.dat")
 
G=   6.67430*pow(10.0,-11.0)
Msun=1.989*pow(10.0,30.0)
Rsun =6.9634*pow(10.0,8.0)
nn=np.zeros((58375,8))
good=0
k=0
for i in range(nr):  
    mass, age, metal1, metal2, logg, teff,Mv,radius= par[i,12], par[i,11], par[i,15], par[i,16],  par[i,9], par[i,8], par[i,5], par[i,14]
    if(mass>0.65 and mass<0.85 and radius>0.1):  
        good+=1
        nn[k,:]=np.array([mass, age, metal1, metal2, logg, teff,Mv,radius])
        #GG=pow(10.0,logg)*0.01
        #gra= G*mass*Msun*(radius*Rsun)**-2.0
        #print ("Gravity:  ",  GG,    gra  )
        fij=open("./Besancon75.dat","a+")
        np.savetxt(fij, nn[k,:].reshape((1,8)),fmt= "%.3f   %.4f    %.3f     %.3f   %.3f    %.1f    %.2f    %.4f")
        k+=1
        fij.close()
    #if(mass>0.65 and mass<0.85 and radius<0.1):  
    #    print ("Error mass,   radius, CL, type:  ",   mass, radius, par[i,6],   par[i,7])   
    #    input("Enter a number ") 
print (good)        


plt.cla()
plt.clf()
plt.figure(figsize=(8, 6))
plt.scatter(nn[:k,0], nn[:k,6],marker='^',c='blue', s=9.0, alpha=1.0)
#plt.plot(tim, Flux,  "r--", lw=1.2)
plt.ylabel(r"$radius$", fontsize=18)
plt.xlabel(r"$mass$", fontsize=18)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.grid("True")
plt.grid(linestyle='dashed')
fig3=plt.gcf()
fig3.savefig("massLumi.jpg", dpi=200)



'''
f0=open("./Besancon75.dat","r")
num0=sum(1 for line in f0)
par=np.zeros(( num0, 8))
par=np.loadtxt("./Besancon75.dat") 
file3=open("Besancon75_sort.dat","w")
file3.close()
massi=np.zeros((num0)); csort=np.zeros((num0))
massi=par[:,0] 
csort=np.argsort(massi); 
for i in range(num0):
    file3=open("Besancon75_sort.dat","a")
    np.savetxt(file3,par[int(csort[i]),:].reshape(1,8),fmt='%.3f   %.4f    %.3f     %.3f   %.3f    %.1f    %.2f    %.4f') 
    file3.close()

################################################################################
par=np.zeros(( num0, 8))
par=np.loadtxt("./Besancon75_sort.dat")     
file4=open("Besancon75_sort2.dat","w")
file4.close()
massr=np.arange(0.65, 0.85, 0.002, dtype=float)

count= np.zeros((len(massr)))   
mass0= np.zeros((len(massr))) 
kk=0
print (massr,  len(massr))
arry=np.zeros((3000, 8))
for i in range(len(massr)):
    k=0
    for j in range(num0): 
        if(abs(par[j,0]-massr[i])<0.0008): 
            arry[k,:]=par[j,:]
            k+=1
            kk+=1
            mass0[i]=par[j,0]
    print ("mass:  ",  massr[i],   k)        
    
    #if(i==0):  count[i]=0
    #else:      
    count[i]=kk # + count[i-1]

    
    
    csort=np.zeros((k))        
    csort=np.argsort(arry[:k,1]) 
    for l in range(k):  
        if(l>0 and arry[int(csort[l]),1]!=arry[int(csort[l-1]),1]):
            file4=open("Besancon75_sort2.dat","a+")
            np.savetxt(file4,arry[int(csort[l]),:].reshape(1,8),fmt='%.3f   %.4f    %.3f     %.3f   %.3f    %.1f    %.2f    %.4f') 
            file4.close()
        else:  count[i]-=1    
################################################################################    
    
print("Count:  ",  count)
print (len(count))    
file4=open("Besancon75_sort2.dat","a+")
np.savetxt(file4,count.reshape(1, len(massr)),fmt='%d, ')
np.savetxt(file4,mass0.reshape(1, len(massr)),fmt='%.3f, ') 
file4.close()    


'''




