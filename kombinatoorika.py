from math import *
import numpy as np  
import pandas as pd
import matplotlib.pyplot as plt 
from scipy.interpolate import spline
from pylab import *
from scipy.integrate import simps

def sellmeier(B1,B2,B3,C1,C2,C3):
    formula = "((" + str(B1) + "*x**2)/(x**2-" + str(C1) + ")+(" + str(B2) + "*x**2)/(x**2-" + str(C2) + ")+(" + str(B3) + "*x**2)/(x**2-" + str(C3) + ")+1)**1/2"

    return formula

#Sellmeieri valemi koefitsendid 
quartz_B1=0.6961663
quartz_B2=0.4079426
quartz_B3=0.8974794
quartz_C1=0.0684043**2
quartz_C2=0.1162414**2
quartz_C3=9.896161**2

BK7_B1=1.03961212
BK7_B2=0.231792344
BK7_B3=1.01046945
BK7_C1=0.00600069867
BK7_C2=0.0200179144
BK7_C3=103.560653

sapphire_B1=1.4313493
sapphire_B2=0.65054713
sapphire_B3=5.3414021
sapphire_C1=0.0726631**2
sapphire_C2=0.1193242**2
sapphire_C3=18.028251**2

CaF2_B1=0.5675888
CaF2_B2=0.4710914
CaF2_B3=3.8484723
CaF2_C1=0.050263605**2
CaF2_C2=0.1003909**2
CaF2_C3=34.649040**2

MgF2_B1=0.48755108
MgF2_B2=0.39875031
MgF2_B3=2.3120353
MgF2_C1=0.04338408**2
MgF2_C2=0.09461442**2
MgF2_C3=23.793604**2

ZnSe_B1=4.45813734
ZnSe_B2=0.467216334
ZnSe_B3=2.89566290
ZnSe_C1=0.200859853**2
ZnSe_C2=0.391371166**2
ZnSe_C3=47.1362108**2




def D_koef_air():
    # Lainepikkus mikromeetrites
    lam = linspace(0.6, 1.6, 1000)
    #B1=0.05792105
    #C1=238.0185
    #B2=0.00167917
    #C2=57.362
    
    #n = 1 + B1/(C1-lam**(-2)) + B2/(C2-lam**(-2))
    
    n = linspace(1,1,1000)
    c = 299792458 * 1e-15  # km/ps
    

    return lam*1e3, -lam*1e3/c * gradient(gradient(n, lam*1e3, edge_order=2), lam*1e3, edge_order=2)


#Sellmeieri valem
def sellmeier(x,B1,B2,B3,C1,C2,C3):
    formula=((B1*x**2)/(x**2-C1)+(B2*x**2)/(x**2-C2)+(B3*x**2)/(x**2-C3)+1)**1/2
    return formula

#Dispersioonikoefitsendi graafik xls-faili ehk andmepunktide põhjal
def disp_xls(filename,title):
    df = pd.read_excel(filename)
    
    wavelength = df.iloc[:,0:1]
    dispersion = df.iloc[:,1:2]
    
    x = []
    y = []
      
    for i in range (0,len(df.index)):
        x.append(df.iloc[i][0])
        
    for i in range (0,len(df.index)):
        y.append(df.iloc[i][1])

        
    plt.plot(x, y, color="#45ad43") 
    plt.title('Fiibri '+filename.split("_",1)[0]+' dispersioon ' + title,fontsize=15)
    
    plt.xlabel('lainepikkus (nm)',fontsize=15)
    plt.xticks(np.arange(min(x), max(x)+1, 100))
    
    plt.ylabel('dispersioon' +'\n'+ '(ps/nm/km)',fontsize=15,rotation=0,labelpad=50)
    plt.show()



#Dispersioonikoefitsendi valem
def D_koef(B1,B2,B3,C1,C2,C3):
    # Lainepikkus mikromeetrites
    lam = linspace(0.6, 1.6, 1000)
    
    n = np.sqrt( 1 + B1 * lam**2 / ( lam**2 - C1 )                       + B2 * lam**2 / ( lam**2 - C2 )                       + B3* lam**2 / ( lam**2 - C3 ))
    
    c = 299792458 * 1e-15  # km/ps

    return lam*1e3, -lam*1e3/c * gradient(gradient(n, lam*1e3, edge_order=2), lam*1e3, edge_order=2)

#Dispersioonikoefitsendi andmed alates valitud lainepikkusest
def D_koef_f(filename):
    df = pd.read_excel(filename)
    
    wavelength = df.iloc[:,0:1]
    dispersion = df.iloc[:,1:2]
    
    x = []
    y = []
      
    for i in range (0,len(df.index)):
        x.append(df.iloc[i][0])
        
    for i in range (0,len(df.index)):
        y.append(df.iloc[i][1])
        
    #Valitud lainepikkus 
    ind = x.index(600)

    return x[ind:],y[ind:]



D_630HP = D_koef_f('630HP_disp.xls')


#muudame kõik listi liikmed täisarvudeks
def inte(lst):
    lst=list(lst)
    for i in range(len(lst)):
        lst[i]=int(lst[i])
    return lst

#loome andmestikule n andmepunkti juurde
def incr_i(data,n):
    x = data[0]
    y = data[1]

    xnew = np.linspace(min(x),max(x),n)
    ynew = spline(x,y,xnew)
    
    return xnew,ynew

#-----------------------------------------------------------------------------------------------------------------------------
#Optimeerimine
Df=D_630HP #<-- uuritav fiiber (mõõteõlal)
lf=35 # <-- maksimum kombinatsiooni pikkus

#potentsiaalsete võrdlusõla klaaside dispersioonikoefitsendid
D_630HP = D_koef_f('630HP_disp.xls')
D_FS=D_koef(quartz_B1,quartz_B2,quartz_B3,quartz_C1,quartz_C2,quartz_C3)
D_BK7=D_koef(BK7_B1,BK7_B2,BK7_B3,BK7_C1,BK7_C2,BK7_C3)
D_sap=D_koef(sapphire_B1,sapphire_B2,sapphire_B3,sapphire_C1,sapphire_C2,sapphire_C3)
D_CaF2=D_koef(CaF2_B1,CaF2_B2,CaF2_B3,CaF2_C1,CaF2_C2,CaF2_C3)
D_MgF2=D_koef(MgF2_B1,MgF2_B2,MgF2_B3,MgF2_C1,MgF2_C2,MgF2_C3)
D_ZnSe=D_koef(ZnSe_B1,ZnSe_B2,ZnSe_B3,ZnSe_C1,ZnSe_C2,ZnSe_C3)#-------------------------6
D_air=D_koef_air()#---------------------------------------------------------------------7

#wl_0 ehk lambda_0 - optimeerime lainepikkuse 633 nanomeetri jaoks
wl_0=685

#Leiame lainepikkusele lambda_0 vastavad dispersioonikoefitsendi väärtused
D_630HP = incr_i(D_630HP,1000)
D_630HP_0=D_630HP[1][inte(D_630HP[0]).index(wl_0)]/1e6
D_FS_0=D_FS[1][inte(D_FS[0]).index(wl_0)]/1e6
D_BK7_0=D_BK7[1][inte(D_BK7[0]).index(wl_0)]/1e6
D_sap_0=D_sap[1][inte(D_sap[0]).index(wl_0)]/1e6
D_CaF2_0=D_CaF2[1][inte(D_CaF2[0]).index(wl_0)]/1e6
D_MgF2_0=D_MgF2[1][inte(D_MgF2[0]).index(wl_0)]/1e6
D_ZnSe_0=D_ZnSe[1][inte(D_ZnSe[0]).index(wl_0)]/1e6
D_air_0=D_air[1][inte(D_air[0]).index(wl_0)]/1e6

n=20 #samm
D_0=[[n*D_FS_0,n*D_BK7_0,n*D_sap_0,n*D_CaF2_0,n*D_MgF2_0,n*D_ZnSe_0,n*D_air_0],['D_FS','D_BK7','D_sap','D_CaF2','D_MgF2','D_ZnSe','D_air']]
#D_0=[[n*D_FS_0,n*D_air_0],['D_FS','D_air']]
print(D_0)
print(D_630HP_0)

import itertools as it
import time
start_time = time.time()

#Klaaside dispersioonikoefitsendi väärtused kombineeritult
D_com=list(it.combinations_with_replacement(D_0[0], lf))
#Vastavalt klaaside nimed
D_com_v=list(it.combinations_with_replacement(D_0[1], lf))

print(len(D_com))
print("kombinatsioonide arv")


#LAMBDA-FUNKTSIOONIGA
#Liidame dispersioonikoefitsendid
D_com_sum=list(map(sum, D_com))

#Leiame parima klaaside kombinatsiooni
#best_fit[0]=dispersioonikoefitsent
#best_fit[1]=kasutatud klaaside nimed (üks kirje iga 1mm kohta)
#best_fit[2]=klaasidega saavutatava dispersioonikoefitsendi erinevus lainepikkusel 633nm fiibri disp.koefitsendist samal lainepikkusel
best_fit=[0,'nimi',50000]
best_fit[0]=min(np.array(D_com_sum), key=lambda x:abs(x-D_630HP_0*419))
best_fit[1]=D_com_v[D_com_sum.index(best_fit[0])]
best_fit[2]=abs(best_fit[0]-D_630HP_0*lf)
print(best_fit)

print("--- %s minutes ---" % str((time.time() - start_time)/60))


klaasid={'D_FS':0,'D_BK7':0,'D_sap':0,'D_CaF2':0,'D_MgF2':0,'D_ZnSe':0,'D_air':0}


best_fit_D2=best_fit

# leiame kui pikk peab iga klaas olema
for i in best_fit_D2[1]:
    for j in klaasid.keys():
        if i==j:
            klaasid[j]+=1
    
print(klaasid)
