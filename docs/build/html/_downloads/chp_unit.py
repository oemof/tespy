from tespy import con, cmp, nwk, hlp
from CoolProp.CoolProp import PropsSI
import math
import matplotlib.pyplot as plt

# %% components

# turbine part
vessel_turb=cmp.vessel(label='vessel_turb',mode='man')
turbine_hp=cmp.turbine(label='turbine_hp',eta_s=0.9)
split=cmp.splitter(label='splitter1')
turbine_lp=cmp.turbine(label='turbine_lp',eta_s=0.9)

# condenser and preheater
condenser=cmp.condenser(label='condenser',dp1=0.95,dp2=0.95,ttd_u=7)
preheater=cmp.condenser(label='preheater',dp1=0.95,dp2=0.99,ttd_u=7)
vessel=cmp.vessel(label='vessel1',mode='man')
merge=cmp.merge(label='merge1')

# feed water
pump=cmp.pump(label='pump',eta_s=0.8,mode='man')
steam_generator=cmp.heat_exchanger_simple(label='steam generator',dp=0.95,mode='man')

# sources and sinks
source=cmp.source(label='source')
sink=cmp.sink(label='sink')

# for cooling water
source_cw=cmp.source(label='source_cw',)
sink_cw=cmp.sink(label='sink_cw')

# %% connections

# turbine part
fs_in=con.connection(source,'out1',vessel_turb,'in1',p=100,T=550)
fs=con.connection(vessel_turb,'out1',turbine_hp,'in1',p=100,m=47,fluid={'water':1})
ext=con.connection(turbine_hp,'out1',split,'in1',p=10)
ext_pre=con.connection(split,'out1',preheater,'in1')
ext_turb=con.connection(split,'out2',turbine_lp,'in1',h0=2800000)

# preheater and condenser
ext_cond=con.connection(preheater,'out1',vessel,'in1',h0=600000)
cond_ws=con.connection(vessel,'out1',merge,'in2')
turb_ws=con.connection(turbine_lp,'out1',merge,'in1')
ws=con.connection(merge,'out1',condenser,'in1')

# feed water
cond=con.connection(condenser,'out1',pump,'in1')
fw_c=con.connection(pump,'out1',preheater,'in2',h0=330000)
fw_w=con.connection(preheater,'out2',steam_generator,'in1')
out=con.connection(steam_generator,'out1',sink,'in1',p=con.ref(fs_in,1,0),h=con.ref(fs_in,1,0))

# cooling water
cw_in=con.connection(source_cw,'out1',condenser,'in2',T=60,p=10,fluid={'water':1})
cw_out=con.connection(condenser,'out2',sink_cw,'in1',T=110)

# %% network

fluids=['water']

nw=nwk.network(fluids=fluids,p='bar',T='C')
nw.add_conns(fs_in,fs,ext,ext_pre,ext_turb,ext_cond,cond_ws,turb_ws,ws,cond,fw_c,fw_w,cw_in,cw_out,out)

# %% busses

power_bus=con.bus('power')
power_bus.add_comp([turbine_hp,-1],[turbine_lp,-1],[pump,-1])

heat_bus=con.bus('heat')
heat_bus.add_comp([condenser,-1])

nw.add_busses(power_bus,heat_bus)


# %% solving

mode='design'

nw.solve(init_file=None,mode=mode)
nw.process_components(mode='post')
nw.save('right_pre_'+mode)

file='right_pre_'+mode+'_results.csv'
mode='offdesign'

fs.set_attr(p=math.nan)
ext.set_attr(p=math.nan)

# overload to partload
m=[50,45,40,35,30]

# temperature level for heating system
t_vl=[80,90,100,110,120]

P={}
Q={}

# iterate over temperatures
for i in t_vl:
    cw_out.set_attr(T=i)
    P[i]=[]
    Q[i]=[]
	# iterate over mass flow
    for j in m:
        fs.set_attr(m=j)

        nw.solve(init_file=file,design_file=file,mode=mode)
        nw.process_components(mode='post')
        P[i]+=[power_bus.P]
        Q[i]+=[heat_bus.P]

# plotting
colors=['#00395b','#74adc1','#bfbfbf','#b54036','#ec6707']

fig, ax = plt.subplots()
j=0
for i in t_vl:
    plt.plot(Q[i],P[i],'.-',Color=colors[j],label='$T_{VL}$ = '+str(i)+' °C', markersize=15,linewidth=2)
    j+=1
ax.set_ylabel('$P$ in W')
ax.set_xlabel('$\dot{Q}$ in W')
plt.legend(loc='upper left')

ax.set_ylim([0, 1e8])
ax.set_xlim([0, 1e8])

plt.show()

fig.savefig('PQ_diagram.svg')