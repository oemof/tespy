from tespy.networks import Network
from tespy.components import (HeatExchanger, Compressor, CycleCloser, Valve, Source, Sink, Merge, DropletSeparator)
from tespy.connections import Connection, Bus
from CoolProp.CoolProp import PropsSI as CPSI
from tespy.tools import ExergyAnalysis
#import plotly.graph_objects as go

wf = 'R1233ZD(E)'           # REFPROP::
si = 'H2O'                  # REFPROP::

# Definition des Netwerks
nw = Network(fluids=[wf, si], T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s', iterinfo=False)

# Definition der Komponenten
GK = HeatExchanger('Gaskühler')
VD = HeatExchanger('Verdampfer')
DR1 = Valve('Drossel 1')
DR2 = Valve('Drossel 2')
KP1 = Compressor('Kompressor 1')
KP2 = Compressor('Kompressor 2')
IWUE = HeatExchanger("Interner Wärmeübertrager")
PH = DropletSeparator('Phasentrenner')
ZU = Merge('Zusammenführung', num_in=2)

# Definition der Quelle, Senke und des Kreislaufzusammenschlusses
se_ein = Source('Senke ein')
se_aus = Sink('Senke aus')

qu_ein = Source('Quelle ein')
qu_aus = Sink('Quelle aus')

KR = CycleCloser('Kreislaufzusammenschluss')

# Verbindungen des Kreislaufs
c21 = Connection(KR, 'out1', GK, 'in1', label="21")
c22 = Connection(GK, 'out1', IWUE, 'in1', label="22")
c23 = Connection(IWUE, 'out1', DR1, 'in1', label="23")
c24 = Connection(DR1, 'out1', PH, 'in1', label="24")
c25 = Connection(PH, 'out1', DR2, 'in1', label="25")
c26 = Connection(DR2, 'out1', VD, 'in2', label="26")
c27 = Connection(VD, 'out2', IWUE, 'in2', label="27")
c28 = Connection(IWUE, 'out2', KP1, 'in1', label="28")
c29 = Connection(KP1, 'out1', ZU, 'in1', label="29")
c30 = Connection(PH, 'out2', ZU, 'in2', label="30")
c31 = Connection(ZU, 'out1', KP2, 'in1', label="31")
c21_cc = Connection(KP2, 'out1', KR, 'in1', label="21_cc")

# Verbindungen der Quelle
c11 = Connection(qu_ein, 'out1', VD, 'in1', label="11")
c12 = Connection(VD, 'out1', qu_aus, 'in1', label="12")

# Verbindungen der Senke
c13 = Connection(se_ein, 'out1', GK, 'in2', label="13")
c14 = Connection(GK, 'out2', se_aus, 'in1', label="14")

nw.add_conns(c21, c22, c23, c24, c25, c26, c27, c28, c29, c30, c31, c21_cc, c11, c12, c13, c14)

# Setzen der Startparameter für die Komponenten
GK.set_attr(pr1=1, pr2=1, Q=-1e7)
IWUE.set_attr(pr1=1, pr2=1)
VD.set_attr(pr1=1, pr2=1)
KP1.set_attr(eta_s=0.76)
KP2.set_attr(eta_s=0.76)

# Setzen Startparameter der Verbindungen des Kreislaufs
h_c22 = CPSI("H", "P", 38 * 1e5, "T", 273.15+165, wf) * 1e-3
c22.set_attr(h=h_c22, p=38)

h_c27 = CPSI("H", "P", 8.334 * 1e5, "T", 273.15+90.1, wf) * 1e-3
c27.set_attr(h=h_c27)

h_c28 = CPSI("H", "P", 8.334 * 1e5, "T", 273.15+155, wf) * 1e-3
c28.set_attr(h=h_c28, p=8.334, fluid={'R1233ZD(E)': 1, 'H2O': 0})

c29.set_attr(p=14.5)

# Setzen Startparameter der Verbindungen der Quelle
c11.set_attr(T=95, p=5, fluid={'R1233ZD(E)': 0, 'H2O': 1})
c12.set_attr(T=90)

# Setzen Startparameter der Verbindungen der Senke
c13.set_attr(T=160, p=20, fluid={'R1233ZD(E)': 0, 'H2O': 1})
c14.set_attr(T=190)

#Lösen des Netzwerks
nw.solve(mode='design')
nw.print_results()

#Setzen der Betriebsparameter
c22.set_attr(h=None, p=43.52)
GK.set_attr(ttd_l=10)
c27.set_attr(h=None, Td_bp=0.1)
c28.set_attr(p=None, h=None)
VD.set_attr(ttd_l=5)
IWUE.set_attr(ttd_u=15)
c29.set_attr(p=16.52)

# Definition der Energieströme
el = Bus('elektrische Leistung')
el.add_comps(
    {'comp': KP1, 'char': 1, 'base': 'bus'},
    {'comp': KP2, 'char': 1, 'base': 'bus'})

wae_zu = Bus('Wärmequelle')
wae_zu.add_comps(
    {'comp': qu_ein, 'base': 'bus'},
    {'comp': qu_aus})

wae_ab = Bus('Wärmesenke')
wae_ab.add_comps(
    {'comp': se_ein, 'base': 'bus'},
    {'comp': se_aus})

nw.add_busses(el, wae_zu, wae_ab)

#Lösen des Netzwerks
nw.solve(mode='design')
nw.print_results()

#Durchführung der Exergianalyse
p_umg = 1
T_umg = 25

# exergy and exergoeconomic analysis
exe_eco_input = {'Gaskühler_Z': 5, 'Drossel 1_Z': 2, 'Drossel 2_Z': 2, 'Phasentrenner_Z': 4,
                 'Kompressor 1_Z': 4, 'Kompressor 2_Z': 4, 'Verdampfer_Z': 4,
                 'Interner Wärmeübertrager_Z': 2,
                 'Zusammenführung_Z': 1,
                 'Quelle ein_c': 0.02, 'Senke ein_c': 0.01, 'elektrische Leistung_c': 0.1}
ean = ExergyAnalysis(nw, E_P=[wae_ab], E_F=[el, wae_zu])
ean.analyse(pamb=p_umg, Tamb=T_umg)
ean.evaluate_exergoeconomics(Exe_Eco_Costs=exe_eco_input, Tamb=T_umg)
ean.print_results(Exe_Eco_An=True)

"""
# Erstellung des Grassmanndiagramms
links, nodes = ean.generate_plotly_sankey_input()
fig = go.Figure(go.Sankey(
    arrangement="snap",
    node={
        "label": nodes,
        'pad': 11,
        'color': 'orange'},
    link=links),
    layout=go.Layout({'width': 1450})
    )
fig.update_layout(
    font_size=20
)
fig.show()
"""
