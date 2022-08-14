import numpy as np
import pygmo as pg

from tespy.components.basics.cycle_closer import CycleCloser
from tespy.components.heat_exchangers.condenser import Condenser
from tespy.components.heat_exchangers.simple import HeatExchangerSimple
from tespy.components.nodes.merge import Merge
from tespy.components.nodes.splitter import Splitter
from tespy.components.piping.valve import Valve
from tespy.components.turbomachinery.pump import Pump
from tespy.components.turbomachinery.turbine import Turbine
from tespy.connections.bus import Bus
from tespy.connections.connection import Connection
from tespy.networks.network import Network
from tespy.tools.optimization import OptimizationProblem


class SamplePlant:
    """Class template for TESPy model usage in optimization module."""
    def __init__(self):

        self.nw = Network(fluids=['water'])
        self.nw.set_attr(p_unit="bar", T_unit="C", h_unit="kJ / kg", iterinfo=False)

        # main cycle components cycle closer
        steam_generator = HeatExchangerSimple("steam generator")
        close_cycle = CycleCloser("cycle closer")

        turbine_hp = Turbine("turbine high pressure")
        turbine_lp = Turbine("turbine low pressure")
        extraction = Splitter("steam extraction splitter", num_out=2)
        preheater = Condenser("feed water preheater")
        valve = Valve("preheater condensate valve")
        waste_steam_merge = Merge("waste steam merge")

        condenser = HeatExchangerSimple("main condenser")
        feed_pump = Pump("feed water pump")

        # Connections

        # main cycle
        c0 = Connection(steam_generator, "out1", close_cycle, "in1", label="0")
        c1 = Connection(close_cycle, "out1", turbine_hp, "in1", label="1")
        c2 = Connection(turbine_hp, "out1", extraction, "in1", label="2")
        c3 = Connection(extraction, "out1", turbine_lp, "in1", label="3")
        c4 = Connection(turbine_lp, "out1", waste_steam_merge, "in1", label="4")
        c5 = Connection(waste_steam_merge, "out1", condenser, "in1", label="5")
        c6 = Connection(condenser, "out1", feed_pump, "in1", label="6")
        c7 = Connection(feed_pump, "out1", preheater, "in2", label="7")
        c8 = Connection(preheater, "out2", steam_generator, "in1", label="8")

        # steam extraction
        c11 = Connection(extraction, "out2", preheater, "in1", label="11")
        c12 = Connection(preheater, "out1", valve, "in1", label="12")
        c13 = Connection(valve, "out1", waste_steam_merge, "in2", label="13")

        self.nw.add_conns(c0, c1, c2, c3, c4, c5, c6, c7, c8, c11, c12, c13)

        # component specifications
        steam_generator.set_attr(pr=0.92)
        turbine_hp.set_attr(eta_s=0.9)
        turbine_lp.set_attr(eta_s=0.9)
        condenser.set_attr(pr=1)
        feed_pump.set_attr(eta_s=0.75)
        preheater.set_attr(ttd_u=5, pr1=1, pr2=0.98)

        # connection specifications
        c1.set_attr(fluid={'water': 1}, p=100, T=600, m=10)
        # pressure at connection 2 will be the parameter to optimize
        c2.set_attr(p=10)
        c6.set_attr(x=0, p=0.1)
        # set liquid state to provide good starting values
        c8.set_attr(state='l')

        power_bus = Bus('power output')
        power_bus.add_comps(
            {'comp': turbine_hp, 'char': 0.97},
            {'comp': turbine_lp, 'char': 0.97},
            {'comp': feed_pump, 'char': 0.97, 'base': 'bus'}
        )
        heat_bus = Bus('heat input')
        heat_bus.add_comps({'comp': steam_generator})
        self.nw.add_busses(power_bus, heat_bus)

        self.nw.solve("design")
        self.stable = "_stable"
        self.nw.save(self.stable)

    def get_param(self, obj, label, parameter):
        """Get the value of a parameter in the network's unit system.

        Parameters
        ----------
        obj : str
            Object to get parameter for (Components/Connections).

        label : str
            Label of the object in the TESPy model.

        parameter : str
            Name of the parameter of the object.

        Returns
        -------
        value : float
            Value of the parameter.
        """
        if obj == 'Components':
            return self.nw.get_comp(label).get_attr(parameter).val
        elif obj == 'Connections':
            return self.nw.get_conn(label).get_attr(parameter).val

    def set_params(self, **kwargs):

        if "Connections" in kwargs:
            for c, params in kwargs["Connections"].items():
                self.nw.get_conn(c).set_attr(**params)

        if "Components" in kwargs:
            for c, params in kwargs["Components"].items():
                self.nw.get_comp(c).set_attr(**params)

    def solve_model(self, **kwargs):
        """
        Solve the TESPy model given the
        """

        self.set_params(**kwargs)

        self.solved = False
        try:
            self.nw.solve("design")
            if self.nw.res[-1] >= 1e-3 or self.nw.lin_dep:
                self.nw.solve("design", init_only=True, init_path=self.stable)
            else:
                # might need more checks here!
                if any(self.nw.results['Condenser']['Q'] > 0):
                    self.solved = False
                else:
                    self.solved = True
        except ValueError as e:
            self.nw.lin_dep = True
            self.nw.solve("design", init_only=True, init_path=self.stable)

    def get_objective(self, objective):
        """
        Get the current objective function evaluation.

        Parameters
        ----------
        objective : str
            Name of the objective function.

        Returns
        -------
        objective_value : float
            Evaluation of the objective function.
        """
        if self.solved:
            return -1 / (
                self.nw.busses['power output'].P.val /
                self.nw.busses['heat input'].P.val
            )
        else:
            return np.nan


plant = SamplePlant()
variables = {"Connections": {"2": {"p": {"min": 0.4, "max": 50}}}}
optimize = OptimizationProblem(plant, variables)

num_ind = 10
num_gen = 15

algo = pg.de()
optimize.run(algo, num_ind, num_gen)
