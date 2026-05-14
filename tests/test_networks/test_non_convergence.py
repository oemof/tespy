
from tespy.components import HeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Turbine
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks import Network


def test_R601_converges_with_linesearch():
    orc = Network()

    caso = Source('cooling air source')
    casi = Sink('cooling air sink')

    turb = Turbine('turbine')
    cond = HeatExchanger('condenser')

    c1 = Connection(caso, 'out1', cond, 'in2', 'c1')
    c2 = Connection(cond, 'out2', casi, 'in1', 'c2')

    b1 = Connection(Source("source"), 'out1', turb, 'in1', 'b1')
    b2 = Connection(turb, 'out1', cond, 'in1', 'b2')
    b3 = Connection(cond, 'out1', Sink("sink"), 'in1', 'b3')

    orc.add_conns(c1, c2, b1, b2, b3)

    # boundary conditions
    air = {"Air":1}
    orcfluid = {"R601": 1}

    lsT = 140 + 273.15
    lsx = 1

    coT = 20 + 273.15
    cop = 1e5
    pr_lossless = 1

    eta_t = 0.85
    ttd_l_c = 15
    t_rise = 10

    c1.set_attr(fluid=air, p=cop, T=coT)
    c2.set_attr(T=Ref(c1, 1, t_rise))
    b1.set_attr(T=lsT, x=lsx, fluid=orcfluid, m=1)
    b3.set_attr(x=0)

    turb.set_attr(eta_s = eta_t)
    cond.set_attr(pr1=pr_lossless, pr2=pr_lossless)
    cond.set_attr(ttd_l=ttd_l_c)

    orc.solve('design', line_search=True)
    orc.assert_convergence()
