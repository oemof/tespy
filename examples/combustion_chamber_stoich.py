
from tespy import con, cmp, nwk, hlp

# %% network

# define full fluid list for the network's variable space
fluid_list = ['TESPy::myAir', 'TESPy::myFuel', 'TESPy::myFuel_fg']
# define unit systems and fluid property ranges
nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
                 p_range=[0.001, 10], T_range=[10, 2000])

# %% components

# sinks & sources
amb = cmp.source('ambient')
sf = cmp.source('fuel')
fg = cmp.sink('flue gas outlet')

# combustion chamber

comb = cmp.combustion_chamber_stoich('stoichiometric combustion chamber')

# %% connections

amb_comb = con.connection(amb, 'out1', comb, 'in1')
sf_comb = con.connection(sf, 'out1', comb, 'in2')
comb_fg = con.connection(comb, 'out1', fg, 'in1')

nw.add_conns(sf_comb, amb_comb, comb_fg)

# %% component parameters

comb.set_attr(fuel={'CH4': 0.96, 'CO2': 0.04},
              air={'Ar': 0.0129, 'N2': 0.7553, 'H2O': 0,
                   'CH4': 0, 'CO2': 0.0004, 'O2': 0.2314},
              fuel_alias='myFuel', air_alias='myAir',
              lamb=3, ti=20000)
# %% connection parameters

# air from abient (ambient pressure and temperature), air composition must be
# stated component wise.
amb_comb.set_attr(T=20, p=1,
                  fluid={'TESPy::myAir': 1, 'TESPy::myFuel': 0,
                         'TESPy::myFuel_fg': 0})

# fuel, pressure must not be stated, as pressure is the same at all inlets and
# outlets of the combustion chamber
sf_comb.set_attr(T=25,
                 fluid={'TESPy::myAir': 0, 'TESPy::myFuel': 1,
                        'TESPy::myFuel_fg': 0})

# %% solving

mode = 'design'
nw.solve(mode=mode)
nw.print_results()

comb.set_attr(path='./LUT')

mode = 'design'
nw.solve(mode=mode)
nw.print_results()
