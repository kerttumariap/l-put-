import scipy.optimize
import numpy as np
import pandas as pd
from math import *
import matplotlib.pyplot as plt 
from scipy.interpolate import spline
from scipy.integrate import simps
from pylab import *
from anneal import Annealer
import random
from scipy import integrate

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

BaF10_B1=1.5851495
BaF10_B2=0.143559385
BaF10_B3=1.08521269
BaF10_C1=0.00926681282
BaF10_C2=0.0424489805
BaF10_C3=105.613573

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


#Dispersioonikoefitsendi valem
def D_koef(B1,B2,B3,C1,C2,C3):
    # Lainepikkus mikromeetrites
    lam = linspace(0.6, 1.6, 1000)
    
    n = np.sqrt( 1 + B1 * lam**2 / ( lam**2 - C1 )                       + B2 * lam**2 / ( lam**2 - C2 )                       + B3* lam**2 / ( lam**2 - C3 ))
    
    c = 299792458 * 1e-15  # km/ps

    return lam*1e3, -lam*1e3/c * gradient(gradient(n, lam*1e3, edge_order=2), lam*1e3, edge_order=2)


def D_koef_air():
    # Lainepikkus mikromeetrites
    lam = linspace(0.6, 1.6, 1000)
    B1=0.05792105
    C1=238.0185
    B2=0.00167917
    C2=57.362
    
    n = 1 + B1/(C1-lam**(-2)) + B2/(C2-lam**(-2))
    
    c = 299792458 * 1e-15  # km/ps

    return lam*1e3, -lam*1e3/c * gradient(gradient(n, lam*1e3, edge_order=2), lam*1e3, edge_order=2)


#Dispersioonikoefitsendi andmete failist lugemine alates valitud lainepikkusest
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


def incr_i(data,n):
    x = data[0]
    y = data[1]

    xnew = np.linspace(min(x),max(x),n)
    ynew = spline(x,y,xnew)
    
    return xnew,ynew

def inte(lst):
    return [int(l) for l in lst]


#----------------------------------------------------------------------------------------------------------------------

# Uuritav kiud 
D_630HP = incr_i(D_koef_f('630HP_disp.xls'),1000)

# Valitud klaasitüübid
D_FS=D_koef(quartz_B1,quartz_B2,quartz_B3,quartz_C1,quartz_C2,quartz_C3)#---------------1
D_BK7=D_koef(BK7_B1,BK7_B2,BK7_B3,BK7_C1,BK7_C2,BK7_C3)#--------------------------------2
D_sap=D_koef(sapphire_B1,sapphire_B2,sapphire_B3,sapphire_C1,sapphire_C2,sapphire_C3)#--3
D_CaF2=D_koef(CaF2_B1,CaF2_B2,CaF2_B3,CaF2_C1,CaF2_C2,CaF2_C3)#-------------------------4
D_MgF2=D_koef(MgF2_B1,MgF2_B2,MgF2_B3,MgF2_C1,MgF2_C2,MgF2_C3)#-------------------------5
D_ZnSe=D_koef(ZnSe_B1,ZnSe_B2,ZnSe_B3,ZnSe_C1,ZnSe_C2,ZnSe_C3)#-------------------------6

klaasid=['Kvartsklaas', 'BK7', 'Safiir', 'Kaltsiumfluoriid', 'Magneesiumfluoriid', 'Tsinkseleniid']


# Vaadeldav lainepikkuste vahemik ühe nanomeetrise sammuga
wl=np.linspace(600,770,170)
Di=D_630HP[1][0:len(wl)]

# Kiu pikkus millimeetrites
li=419 # 419 mm 
# Teisendame dispersioonikoefitsendi ühikutes ühiku kilomeeter millimeetriteks, sest vaatleme kiudu, mis on 419mm
Di=np.divide(Di,1e6)

# Valitud klaaside dispersioonikoefitsentide väärtused valitud lainepikkuste vahemikus wl
D=[D_FS[1][0:len(wl)], D_BK7[1][0:len(wl)], D_sap[1][0:len(wl)], D_CaF2[1][0:len(wl)], D_MgF2[1][0:len(wl)],D_ZnSe[1][0:len(wl)]]
D=np.asarray(D)
# Teisendame dispersioonikoefitsendi ühikutes ühiku kilomeeter millimeetriteks
D=np.divide(D,1e6)

# Valime algväärtuse, kus kvartsklaasi pikkus on 424 mm ning teiste valitud klaaside pikkused on 0 mm
state=[424, 0, 0, 0, 0, 0]
state=np.asarray(state)

# Minimeeritav funktsioon
def energy(state):
    y1=np.power(np.multiply(np.multiply(Di,li),1e2),1)
    y2=np.power(np.multiply(np.multiply(D[0],state[0])+np.multiply(D[1],state[1])+np.multiply(D[2],state[2])+np.multiply(D[3],state[3])+np.multiply(D[4],state[4])+np.multiply(D[5],state[5]),1e2),1)
    eq=y1-y2
    eq=np.multiply(np.sum(np.abs(eq)),1)
    return eq


print(energy(state))

# Valime ühe klaasi suurima võimaliku pikkuse (700 mm) ja muudame käesolevat väärtust vastavalt sellele, n-ö kontrollitult juhuslikult
def move(state):
    i = random.randint(0,5)
    if state[i]<700 and state[i]>0: 
        state[i] = random.randint(state[i]-1,state[i]+1)
    else:
        state[i] = random.randint(0,700)
        

import time
start_time = time.time()

annealer = Annealer(energy, move)
schedule = annealer.auto(state, minutes=5)
state, e = annealer.anneal(state, schedule['tmax'], schedule['tmin'], schedule['steps'], updates=1000)
state, e = annealer.anneal(state, 10000, 0.1, 100000, updates=30)

# Parim klaaside kombinatsioon
for i in range(len(state)): 
    print(klaasid[i])
    print(state[i])

# Minimeeritava funktsiooni väärtus, mis vastab parimale klaasi kombinatsioonile
print('Minimeeritava funktsiooni väärtus:')
print(energy(state))


np.savetxt("sub_uus.csv", state, delimiter=",")

print("--- %s minutes ---" % str(float(time.time() - start_time)/60.0))
