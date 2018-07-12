from tespy import nwk, con, cmp, hlp
import numpy as np
from matplotlib import pyplot as plt

nw = nwk.network(['water'], p_unit='bar', T_unit='C', h_unit='kJ / kg',
                            p_range=[0.5, 2], T_range=[10, 100],
                            h_range=[5,500])

# %% components
pipe = cmp.pipe('pipe')
pipe2 = cmp.pipe('pipe2')
pipe3 = cmp.pipe('pipe3')
pipe4 = cmp.pipe('pipe4')
pipe5 = cmp.pipe('pipe5')
pipe6 = cmp.pipe('pipe6')
sink = cmp.sink('sink')
source = cmp.source('source')

# %% connections

a = con.connection(source, 'out1', pipe, 'in1')
b = con.connection(pipe, 'out1', pipe2, 'in1')
c = con.connection(pipe2, 'out1', pipe3, 'in1')
d = con.connection(pipe3, 'out1', pipe4, 'in1')
e = con.connection(pipe4, 'out1', pipe5, 'in1')
f = con.connection(pipe5, 'out1', pipe6, 'in1')
g = con.connection(pipe6, 'out1', sink, 'in1')

nw.add_conns(a, b, c, d, e, f, g)

# %% connection parameters

a.set_attr(h=40, fluid={'water': 1}, p=1, m=10)
b.set_attr(h=40)
c.set_attr(h=40)
d.set_attr(h=40)
e.set_attr(h=40)
f.set_attr(h=40)
g.set_attr(h=40)


# %% component parameters

pipe.set_attr(pr=0.99, ks=1e-4, L=20, D='var')
pipe2.set_attr(pr=0.95, ks=1e-4, L=20, D='var')
pipe3.set_attr(pr=0.95, ks=1e-4, L=20, D='var')

pipe4.set_attr(pr=0.99, ks=1e-4, L=30, D='var')
pipe5.set_attr(pr=0.95, ks=1e-4, L=20, D='var')
pipe6.set_attr(pr=0.95, ks=1e-4, L=20, D='var')
pipe.D.val=10
pipe2.D.val=10
pipe3.D.val=10

pipe4.D.val=10
pipe5.D.val=10
pipe6.D.val=10

# %% solve

nw.solve(mode='design')
nw.print_results()
nw.save('test')
