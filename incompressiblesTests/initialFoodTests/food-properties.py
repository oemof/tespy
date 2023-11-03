# %%
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 22:41:28 2022

@author: mrk
"""

import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal

print("loaded modules")

# %%

Temps = np.linspace(-40,150,50)

cp = {'Protein' : [],
      'Fat' : [],
      'Carbohydrate' : [],
      'Fiber' : [],
      'Ash' : [],
      'Water' : []
      }
k = {'Protein' : [],
      'Fat' : [],
      'Carbohydrate' : [],
      'Fiber' : [],
      'Ash' : [],
      'Water' : []
      }
a = {'Protein' : [],
      'Fat' : [],
      'Carbohydrate' : [],
      'Fiber' : [],
      'Ash' : [],
      'Water' : []
      }
d = {'Protein' : [],
      'Fat' : [],
      'Carbohydrate' : [],
      'Fiber' : [],
      'Ash' : [],
      'Water' : []
      }
a2 = {'Protein' : [],
      'Fat' : [],
      'Carbohydrate' : [],
      'Fiber' : [],
      'Ash' : [],
      'Water' : []
      }

#Specific heat, kJ/(kg·K) 
for T in Temps:
    cp['Protein']      += [2.0082 + 1.2089 * 1e-3*T - 1.3129 * 1e-6*T**2]
    cp['Fat']          += [1.9842 + 1.4733 * 1e-3*T - 4.8008 * 1e-6*T**2]
    cp['Carbohydrate'] += [1.5488 + 1.9625 * 1e-3*T - 5.9399 * 1e-6*T**2]
    cp['Fiber']        += [1.8459 + 1.8306 * 1e-3*T - 4.6509 * 1e-6*T**2]
    cp['Ash']          += [1.0926 + 1.8896 * 1e-3*T - 3.6817 * 1e-6*T**2]    
    cp['Water']        += [4.1289 - 9.0864 * 1e-5*T + 5.4731 * 1e-6*T**2]

    k['Protein']       += [1.7881 * 1e-1 + 1.1958 * 1e-3*T - 2.7178 * 1e-6*T**2]
    k['Fat']           += [1.8071 * 1e-1 - 2.7604 * 1e-4*T - 1.7749 * 1e-7*T**2]
    k['Carbohydrate']  += [2.0141 * 1e-1 + 1.3874 * 1e-3*T - 4.3312 * 1e-6*T**2]
    k['Fiber']         += [1.8331 * 1e-1 + 1.2497 * 1e-3*T - 3.1683 * 1e-6*T**2]
    k['Ash']           += [3.2962 * 1e-1 + 1.4011 * 1e-3*T - 2.9069 * 1e-6*T**2]
    k['Water']         += [5.7109 * 1e-1 + 1.7625 * 1e-3*T - 6.7036 * 1e-6*T**2]
        
    a['Protein']       += [6.8714 * 1e-8 + 4.7578 * 1e-10*T - 1.4646 * 1e-12*T**2]
    a['Fat']           += [9.8777 * 1e-8 - 1.2569 * 1e-11*T - 3.8286 * 1e-14*T**2]
    a['Carbohydrate']  += [8.0842 * 1e-8 + 5.3052 * 1e-10*T - 2.3218 * 1e-12*T**2]
    a['Fiber']         += [7.3976 * 1e-8 + 5.1902 * 1e-10*T - 2.2202 * 1e-12*T**2]
    a['Ash']           += [1.2461 * 1e-7 + 3.7321 * 1e-10*T - 1.2244 * 1e-12*T**2]
    a['Water']         += [1.3168 * 1e-7 + 6.2477 * 1e-10*T - 2.4022 * 1e-12*T**2]
    
    d['Protein']       += [1.3299 * 1e3 - 5.1840 * 1e-1*T]
    d['Fat']           += [9.2559 * 1e2 - 4.1757 * 1e-1*T]
    d['Carbohydrate']  += [1.5991 * 1e3 - 3.1046 * 1e-1*T]
    d['Fiber']         += [1.3115 * 1e3 - 3.6589 * 1e-1*T]
    d['Ash']           += [2.4238 * 1e3 - 2.8063 * 1e-1*T]
    d['Water']         += [9.9718 * 1e2 + 3.1439 * 1e-3*T - 3.7574 * 1e-3*T**2]
    
a2['Protein']       = [k/(cp*d) for k,cp,d in zip(k['Protein'],cp['Protein'],d['Protein'])]
a2['Fat']           = [k/(cp*d) for k,cp,d in zip(k['Fat'],cp['Fat'],d['Fat'])]
a2['Carbohydrate']  = [k/(cp*d) for k,cp,d in zip(k['Carbohydrate'],cp['Carbohydrate'],d['Carbohydrate'])]
a2['Fiber']         = [k/(cp*d) for k,cp,d in zip(k['Fiber'],cp['Fiber'],d['Fiber'])]
a2['Ash']           = [k/(cp*d) for k,cp,d in zip(k['Ash'],cp['Ash'],d['Ash'])]
a2['Water']         = [k/(cp*d) for k,cp,d in zip(k['Water'],cp['Water'],d['Water'])]
    
    
fig, ax = plt.subplots(2, 2, figsize=(16, 8))
ax = ax.flatten()
[a.grid() for a in ax]
[a.set_xlabel('temperature, C') for a in ax]

for i,key in enumerate(cp):
    ax[0].scatter(Temps, cp[key],label=str(key))
    
for i,key in enumerate(k):
    ax[1].scatter(Temps, k[key],label=str(key))

for i,key in enumerate(a):
    ax[2].scatter(Temps, [a/a2*1000 for a,a2 in zip(a[key],a2[key])],label=str(key))

for i,key in enumerate(d):
    ax[3].scatter(Temps, d[key],label=str(key))

[a.legend() for a in ax]
ax[0].set_ylabel('cp')
ax[1].set_ylabel('k')
ax[2].set_ylabel('a')
ax[3].set_ylabel('d')

fig.show


np.set_printoptions(precision=5,suppress=False)

Temps
print('\n')


for key in [cp,k,d]:
    print(np.array(key['Protein']))
print('\n')
for key in [cp,k,d]:    
    print(np.array(key['Fat']))
print('\n')
for key in [cp,k,d]:    
    print(np.array(key['Carbohydrate']))
print('\n')
for key in [cp,k,d]:    
    print(np.array(key['Fiber']))
print('\n')
for key in [cp,k,d]:    
    print(np.array(key['Ash']))
print('\n')
for key in [cp,k,d]:    
    print(np.array(key['Water']))

# %%


CP.get_global_param_string("version")
fl = CP.get_global_param_string("FluidsList")


CP_cp = {'Protein' : [],
      'Fat' : [],
      'Carbohydrate' : [],
      'Fiber' : [],
      'Ash' : [],
      'Water' : []
      }
CP_k = {'Protein' : [],
      'Fat' : [],
      'Carbohydrate' : [],
      'Fiber' : [],
      'Ash' : [],
      'Water' : []
      }
CP_d = {'Protein' : [],
      'Fat' : [],
      'Carbohydrate' : [],
      'Fiber' : [],
      'Ash' : [],
      'Water' : []
      }

Temps[0] = Temps[0] +1e-4
Temps[-1] = Temps[-1] -1e-4


#Specific heat, kJ/(kg·K) 
for T in Temps:
    CP_cp['Protein']      += [CP.PropsSI('C','T',T+273.15,'P',101325,'INCOMP::FoodProtein')]
    CP_cp['Fat']          += [CP.PropsSI('C','T',T+273.15,'P',101325,'INCOMP::FoodFat')]
    CP_cp['Carbohydrate'] += [CP.PropsSI('C','T',T+273.15,'P',101325,'INCOMP::FoodCarbohydrate')]
    CP_cp['Fiber']        += [CP.PropsSI('C','T',T+273.15,'P',101325,'INCOMP::FoodFiber')]
    CP_cp['Ash']          += [CP.PropsSI('C','T',T+273.15,'P',101325,'INCOMP::FoodAsh')]
    CP_cp['Water']        += [CP.PropsSI('C','T',T+273.15,'P',101325,'INCOMP::FoodWater')]

    CP_k['Protein']       += [CP.PropsSI('L','T',T+273.15,'P',101325,'INCOMP::FoodProtein')]
    CP_k['Fat']           += [CP.PropsSI('L','T',T+273.15,'P',101325,'INCOMP::FoodFat')]
    CP_k['Carbohydrate']  += [CP.PropsSI('L','T',T+273.15,'P',101325,'INCOMP::FoodCarbohydrate')]
    CP_k['Fiber']         += [CP.PropsSI('L','T',T+273.15,'P',101325,'INCOMP::FoodFiber')]
    CP_k['Ash']           += [CP.PropsSI('L','T',T+273.15,'P',101325,'INCOMP::FoodAsh')]
    CP_k['Water']         += [CP.PropsSI('L','T',T+273.15,'P',101325,'INCOMP::FoodWater')]
        
    CP_d['Protein']       += [CP.PropsSI('D','T',T+273.15,'P',101325,'INCOMP::FoodProtein')]
    CP_d['Fat']           += [CP.PropsSI('D','T',T+273.15,'P',101325,'INCOMP::FoodFat')]
    CP_d['Carbohydrate']  += [CP.PropsSI('D','T',T+273.15,'P',101325,'INCOMP::FoodCarbohydrate')]
    CP_d['Fiber']         += [CP.PropsSI('D','T',T+273.15,'P',101325,'INCOMP::FoodFiber')]
    CP_d['Ash']           += [CP.PropsSI('D','T',T+273.15,'P',101325,'INCOMP::FoodAsh')]
    CP_d['Water']         += [CP.PropsSI('D','T',T+273.15,'P',101325,'INCOMP::FoodWater')]
    

   
fig, ax = plt.subplots(2, 2, figsize=(16, 8))
ax = ax.flatten()
[a.grid() for a in ax]
[a.set_xlabel('temperature, C') for a in ax]

for i,key in enumerate(CP_cp):
    ax[0].scatter(Temps, CP_cp[key],label=str(key))
    
for i,key in enumerate(CP_k):
    ax[1].scatter(Temps, CP_k[key],label=str(key))

#    ax[2].scatter(Temps, [a/a2*1000 for a,a2 in zip(a[key],a2[key])],label=str(key))

for i,key in enumerate(CP_d):
    ax[3].scatter(Temps, CP_d[key],label=str(key))

[a.legend() for a in ax]
ax[0].set_ylabel('cp')
ax[1].set_ylabel('k')
ax[3].set_ylabel('d')



plt.show(block=True)



