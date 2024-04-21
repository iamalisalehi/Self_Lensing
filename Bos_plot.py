import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import warnings
warnings.filterwarnings("ignore")
rcParams["font.size"] = 13
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
################################################################################
nm=int(7319)
data=np.zeros((nm, 9))
data=np.loadtxt('Besancon_sort2_Limb.dat')
nam=['mass', 'age', 'metal' , 'logg', 'teff', 'Mv', 'radius', 'limb', 'CL']

ms=0
sg=0
dat1=np.zeros((nm, 9))
dat2=np.zeros((nm, 9))
for i in range(nm):  
    if(data[i,8]==5 and data[i,6]>2.0 and data[i,6]<5.0):  
        dat1[ms,:]= data[i,:]## main-sequence
        ms+=1
    if(data[i,8]!=5 and data[i,6]>2.0 and data[i,6]<5.0):   
        dat2[sg,:]= data[i,:]## ginat
        sg+=1
    
print  ("Main-sequence fraction:  ",  ms*100.0/nm)    
print  ("Subgiant      fraction:  ",  sg*100.0/nm,   ms+sg)        
    


'''
for i in range(9):
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    ax= plt.gca()              
    plt.hist(data[:,i],30,histtype='bar',ec='darkgreen',facecolor='green',alpha=0.7, rwidth=1.5)
    y_vals = ax.get_yticks()
    ax.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/nm)) for x in y_vals]) 
    y_vals = ax.get_yticks()
    plt.ylim([np.min(y_vals), np.max(y_vals)])
    ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=19,labelpad=0.1)
    ax.set_xlabel(str(nam[i]),fontsize=19,labelpad=0.1)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig=plt.gcf()
    fig.savefig("./Histo{0:d}.jpg".format(i),dpi=200)
print ("****  All histos are plotted *****************************")   
'''
################################################################################
#2.044 0.0923 -0.201 3.51 7617.0 0.61 3.488 0.4431 3.0
#2.123 0.1014 0.07 3.46 7519.0 0.43 3.819 0.4115 3.0
#2.127 0.1187 -0.014 3.67 7578.0 0.75 3.65 0.41 3.0

target=np.array([2.123 , 0.1014, 0.07 ,3.46,7519.0,0.43, 3.819, 0.4115, 3.0])

plt.cla()
plt.clf()
plt.figure(figsize=(8, 6))

plt.scatter(dat1[:ms,0], dat1[:ms,1], s=dat1[:ms,6]*dat1[:ms,6]*2.5, marker='s', c=dat1[:ms,7],  alpha=0.7)
plt.scatter(dat2[:sg,0], dat2[:sg,1], s=dat2[:sg,6]*dat2[:sg,6]*2.5, marker='o', c=dat2[:sg,7],  alpha=0.7)
plt.scatter(target[0], target[1], s=target[6]*target[6]*10.5, marker='x', c="m", alpha=0.99 )
cbar=plt.colorbar(pad=0.02)
cbar.set_label(r"$\rm{Limb}-\rm{Darkening}~\rm{coefficient}$", fontsize=18, rotation=270, labelpad=15.0)
cbar.ax.tick_params(labelsize=16)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
#plt.clim(0.38, 0.57)
#plt.yscale('log')
plt.xlim(1.5, 2.37)
plt.ylim(0.0, 1.78)
plt.grid("True")
plt.xlabel(r"$\rm{Mass}(M_{\odot})$", fontsize=18)
plt.ylabel(r"$\rm{Age}(\rm{Gyrs})$", fontsize=18)
plt.grid(linestyle='dashed')
fig3=plt.gcf()
fig3.tight_layout(pad=0.0)
fig3.savefig("./MassAge.jpg", dpi=200)
print(">>>>> The best u, Astar, ...  was made <<<<<<<<<<<<<<<")
'''
file4=open("test.dat","w")
for i in range(nm): 
    if(data[i,1]>1.0): 
        np.savetxt(file4,data[i,:].reshape(1, 8),fmt='%.3f   %.4f    %.3f    %.3f    %.1f    %.2f    %.4f   %.4f')
'''

