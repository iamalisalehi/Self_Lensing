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
################################################################################


fd=open("./Sector6SAPdata.dat","r")
n1=sum(1 for line in fd)
sap1= np.zeros((n1,3))
sap1=np.loadtxt("./Sector6SAPdata.dat") 
plt.cla()
plt.clf()
plt.figure(figsize=(8,3))
plt.errorbar(sap1[:,0],sap1[:,1], yerr=sap1[:,2],fmt=".",markersize=4.2,color='g',ecolor='darkgreen',elinewidth=0.25,capsize=0)
plt.ylabel(r"$\rm{Flux}(e^{-1}/\rm{s})$", fontsize=18)
plt.xlabel(r"$\rm{Time}-2457000[\rm{BTJD}~\rm{days}]$", fontsize=18)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.text(1485, 5690, r"$\rm{Sector}~6$", fontsize=18)
plt.grid("True")
plt.grid(linestyle='dashed')
fig=plt.gcf()
fig.tight_layout(pad=0.2)
fig.savefig("./lightSector6.jpg",dpi=200)
print("Lightcurve_ini is plotted!!!!")

################################################################################

fd=open("./Sector32SAPdata.dat","r")
n1=sum(1 for line in fd)
sap1= np.zeros((n1,3))
sap1=np.loadtxt("./Sector32SAPdata.dat") 
plt.cla()
plt.clf()
plt.figure(figsize=(8,3))
plt.errorbar(sap1[:,0],sap1[:,1], yerr=sap1[:,2],fmt=".",markersize=4.2,color='g',ecolor='darkgreen',elinewidth=0.25,capsize=0)
plt.ylabel(r"$\rm{Flux}(e^{-1}/\rm{s})$", fontsize=18)
plt.xlabel(r"$\rm{Time}-2457000[\rm{BTJD}~\rm{days}]$", fontsize=18)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.text(2193, 6180, r"$\rm{Sector}~32$", fontsize=18)
plt.grid("True")
plt.grid(linestyle='dashed')
fig=plt.gcf()
fig.tight_layout(pad=0.2)
fig.savefig("./lightSector32.jpg",dpi=200)
print("Lightcurve_ini is plotted!!!!")


################################################################################

fd=open("./Sector33SAPdata.dat","r")
n1=sum(1 for line in fd)
sap1= np.zeros((n1,3))
sap1=np.loadtxt("./Sector33SAPdata.dat") 
plt.cla()
plt.clf()

plt.errorbar(sap1[:,0],sap1[:,1], yerr=sap1[:,2],fmt=".",markersize=4.2,color='g',ecolor='darkgreen',elinewidth=0.25,capsize=0)
plt.ylabel(r"$\rm{Flux}(e^{-1}/\rm{s})$", fontsize=18)
plt.xlabel(r"$\rm{Time}-2457000[\rm{BTJD}~\rm{days}]$", fontsize=18)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.text(2222, 6650, r"$\rm{Sector}~33$", fontsize=18)
plt.grid("True")
plt.grid(linestyle='dashed')
fig=plt.gcf()
fig.tight_layout(pad=0.2)
fig.savefig("./lightSector33.jpg",dpi=200)
print("Lightcurve_ini is plotted!!!!")




























