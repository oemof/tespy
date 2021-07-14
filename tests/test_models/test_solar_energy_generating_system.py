# -*- coding: utf-8

"""Module for testing a tespy simulation vs results from a different simulator.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_models/test_solar_energy_generating_system.py

SPDX-License-Identifier: MIT
"""
from tespy.components import Compressor
from tespy.components import Condenser
from tespy.components import CycleCloser
from tespy.components import Drum
from tespy.components import HeatExchanger
from tespy.components import Merge
from tespy.components import ParabolicTrough
from tespy.components import Pump
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import Turbine
from tespy.components import Valve
from tespy.connections import Bus
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks import Network
from tespy.tools import ExergyAnalysis


class TestSEGS:

    def setup(self):
        """
        Full model validation of SEGS model in TESPy vs. EBSILON.

        Find original models at https://github.com/fwitte/SEGS_exergy.
        """
        # specification of ambient state
        self.pamb = 1.013
        self.Tamb = 25

        # setting up network
        self.nw = Network(fluids=['water', 'INCOMP::TVP1', 'air'])
        self.nw.set_attr(
            T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s',
            s_unit="kJ / kgK")

        # components definition
        air_in = Source('Ambient air source', fkt_group='CW')
        air_out = Sink('Ambient air sink', fkt_group='CW')

        closer_pt = CycleCloser('Cycle closer pt', fkt_group='SF')
        pt = ParabolicTrough('Parabolic trough', fkt_group='SF')
        ptpump = Pump('HTF pump', fkt_group='SF')

        closer = CycleCloser('Cycle closer power cycle', fkt_group='SG')

        eco = HeatExchanger('Economizer', fkt_group='SG')
        eva = HeatExchanger('Evaporator', fkt_group='SG')
        sup = HeatExchanger('Superheater', fkt_group='SG')
        drum = Drum('Drum', fkt_group='SG')

        reh = HeatExchanger('Reheater', fkt_group='RH')

        hpt1 = Turbine('HP turbine 1', fkt_group='HPT')
        hpt2 = Turbine('HP turbine 2', fkt_group='HPT')
        lpt1 = Turbine('LP turbine 1', fkt_group='LPT')
        lpt2 = Turbine('LP turbine 2', fkt_group='LPT')
        lpt3 = Turbine('LP turbine 3', fkt_group='LPT')
        lpt4 = Turbine('LP turbine 4', fkt_group='LPT')
        lpt5 = Turbine('LP turbine 5', fkt_group='LPT')

        cond = Condenser('Condenser', fkt_group='CW')
        condpump = Pump('Condenser pump', fkt_group='CW')
        fwt = Merge('Feedwater tank', num_in=3, fkt_group='LPP')
        fwp = Pump('Feedwater pump', fkt_group='FWP')
        cwp = Pump('Cooling water pump', fkt_group='CW')
        closer_cw = CycleCloser('Cycle closer cw', fkt_group='CW')
        ct = HeatExchanger('Cooling tower', fkt_group='CW')
        fan = Compressor('Cooling tower fan', fkt_group='CW')

        sp1 = Splitter('Splitter 1', fkt_group='HPT')
        sp2 = Splitter('Splitter 2', fkt_group='HPT')
        sp3 = Splitter('Splitter 3', fkt_group='LPT')
        sp4 = Splitter('Splitter 4', fkt_group='LPT')
        sp5 = Splitter('Splitter 5', fkt_group='LPT')
        sp6 = Splitter('Splitter 6', fkt_group='LPT')
        sp7 = Splitter('Splitter 7', fkt_group='SF')

        m1 = Merge('Merge 1', fkt_group='CW')
        m2 = Merge('Merge 2', fkt_group='HPP')
        m3 = Merge('Merge 3', fkt_group='LPP')
        m4 = Merge('Merge 4', fkt_group='LPP')
        m5 = Merge('Merge 5', fkt_group='SF')

        v1 = Valve('Valve 1', fkt_group='HPP')
        v2 = Valve('Valve 2', fkt_group='HPP')
        v3 = Valve('Valve 3', fkt_group='LPP')
        v4 = Valve('Valve 4', fkt_group='LPP')
        v5 = Valve('Valve 5', fkt_group='LPP')

        hppre1 = Condenser('High pressure preheater 1', fkt_group='HPP')
        hppre2 = Condenser('High pressure preheater 2', fkt_group='HPP')
        hppre1_sub = HeatExchanger(
            'High pressure preheater 1 subcooling', fkt_group='HPP')
        hppre2_sub = HeatExchanger(
            'High pressure preheater 2 subcooling', fkt_group='HPP')

        lppre1 = Condenser('Low pressure preheater 1', fkt_group='LPP')
        lppre2 = Condenser('Low pressure preheater 2', fkt_group='LPP')
        lppre3 = Condenser('Low pressure preheater 3', fkt_group='LPP')
        lppre1_sub = HeatExchanger(
            'Low pressure preheater 1 subcooling', fkt_group='LPP')
        lppre2_sub = HeatExchanger(
            'Low pressure preheater 2 subcooling', fkt_group='LPP')
        lppre3_sub = HeatExchanger(
            'Low pressure preheater 3 subcooling', fkt_group='LPP')

        # connections definition
        # power cycle
        c1 = Connection(sup, 'out2', closer, 'in1', label='1')
        c2 = Connection(closer, 'out1', hpt1, 'in1', label='2')
        c3 = Connection(hpt1, 'out1', sp1, 'in1', label='3')
        c4 = Connection(sp1, 'out1', hpt2, 'in1', label='4')
        c5 = Connection(hpt2, 'out1', sp2, 'in1', label='5')
        c6 = Connection(sp2, 'out1', reh, 'in2', label='6')
        c7 = Connection(reh, 'out2', lpt1, 'in1', label='7')
        c8 = Connection(lpt1, 'out1', sp3, 'in1', label='8')
        c9 = Connection(sp3, 'out1', lpt2, 'in1', label='9')
        c10 = Connection(lpt2, 'out1', sp4, 'in1', label='10')
        c11 = Connection(sp4, 'out1', lpt3, 'in1', label='11')
        c12 = Connection(lpt3, 'out1', sp5, 'in1', label='12')
        c13 = Connection(sp5, 'out1', lpt4, 'in1', label='13')
        c14 = Connection(lpt4, 'out1', sp6, 'in1', label='14')
        c15 = Connection(sp6, 'out1', lpt5, 'in1', label='15')
        c16 = Connection(lpt5, 'out1', m1, 'in1', label='16')
        c17 = Connection(m1, 'out1', cond, 'in1', label='17')
        c18 = Connection(cond, 'out1', condpump, 'in1', label='18')
        c19 = Connection(condpump, 'out1', lppre1, 'in2', label='19')
        # c19 = Connection(condpump, 'out1', lppre1_sub, 'in2', label='19')
        # c20 = Connection(lppre1_sub, 'out2', lppre1, 'in2', label='20')
        c21 = Connection(lppre1, 'out2', lppre2, 'in2', label='21')
        # c21 = Connection(lppre1, 'out2', lppre2_sub, 'in2', label='21')
        # c22 = Connection(lppre2_sub, 'out2', lppre2, 'in2', label='22')
        c23 = Connection(lppre2, 'out2', lppre3, 'in2', label='23')
        # c23 = Connection(lppre2, 'out2', lppre3_sub, 'in2', label='23')
        # c24 = Connection(lppre3_sub, 'out2', lppre3, 'in2', label='24')
        c25 = Connection(lppre3, 'out2', fwt, 'in1', label='25')
        c26 = Connection(fwt, 'out1', fwp, 'in1', label='26')
        c27 = Connection(fwp, 'out1', hppre1, 'in2', label='27')
        c29 = Connection(hppre1, 'out2', hppre2, 'in2', label='29')
        c31 = Connection(hppre2, 'out2', eco, 'in2', label='31')

        c36 = Connection(sp1, 'out2', hppre2, 'in1', label='36')
        c37 = Connection(hppre2, 'out1', v1, 'in1', label='37')
        c39 = Connection(v1, 'out1', m2, 'in2', label='39')
        c40 = Connection(sp2, 'out2', m2, 'in1', label='40')
        c41 = Connection(m2, 'out1', hppre1, 'in1', label='41')
        c42 = Connection(hppre1, 'out1', v2, 'in1', label='42')
        c44 = Connection(v2, 'out1', fwt, 'in2', label='44')
        c45 = Connection(sp3, 'out2', fwt, 'in3', label='45')
        c46 = Connection(sp4, 'out2', lppre3, 'in1', label='46')
        c47 = Connection(lppre3, 'out1', v3, 'in1', label='47')
        # c47 = Connection(lppre3, 'out1', lppre3_sub, 'in1', label='47')
        # c48 = Connection(lppre3_sub, 'out1', v3, 'in1', label='48')
        c49 = Connection(v3, 'out1', m3, 'in1', label='49')
        c50 = Connection(sp5, 'out2', m3, 'in2', label='50')
        c51 = Connection(m3, 'out1', lppre2, 'in1', label='51')
        c52 = Connection(lppre2, 'out1', v4, 'in1', label='52')
        # c52 = Connection(lppre2, 'out1', lppre2_sub, 'in1', label='52')
        # c53 = Connection(lppre2_sub, 'out1', v4, 'in1', label='53')
        c54 = Connection(v4, 'out1', m4, 'in2', label='54')
        c55 = Connection(sp6, 'out2', m4, 'in1', label='55')
        c56 = Connection(m4, 'out1', lppre1, 'in1', label='56')
        c57 = Connection(lppre1, 'out1', v5, 'in1', label='57')
        # c57 = Connection(lppre1, 'out1', lppre1_sub, 'in1', label='57')
        # c58 = Connection(lppre1_sub, 'out1', v5, 'in1', label='58')
        c59 = Connection(v5, 'out1', m1, 'in2', label='59')

        # components from subsystem
        c32 = Connection(eco, 'out2', drum, 'in1', label='32')
        c33 = Connection(drum, 'out1', eva, 'in2', label='33')
        c34 = Connection(eva, 'out2', drum, 'in2', label='34')
        c35 = Connection(drum, 'out2', sup, 'in2', label='35')
        c73 = Connection(sup, 'out1', eva, 'in1', label='73')
        c74 = Connection(eva, 'out1', eco, 'in1', label='74')

        # cooling water
        c60 = Connection(cond, 'out2', closer_cw, 'in1', label='60')
        c61 = Connection(closer_cw, 'out1', ct, 'in1', label='61')
        c62 = Connection(ct, 'out1', cwp, 'in1', label='62')
        c63 = Connection(cwp, 'out1', cond, 'in2', label='63')

        # cooling tower
        c64 = Connection(air_in, 'out1', fan, 'in1', label='64')
        c65 = Connection(fan, 'out1', ct, 'in2', label='65')
        c66 = Connection(ct, 'out2', air_out, 'in1', label='66')

        # parabolic trough cycle
        c70 = Connection(pt, 'out1', closer_pt, 'in1', label='67')
        c71 = Connection(closer_pt, 'out1', sp7, 'in1', label='71')
        c72 = Connection(sp7, 'out1', sup, 'in1', label='72')
        c75 = Connection(eco, 'out1', m5, 'in1', label='75')
        c76 = Connection(sp7, 'out2', reh, 'in1', label='76')
        c77 = Connection(reh, 'out1', m5, 'in2', label='77')
        c78 = Connection(m5, 'out1', ptpump, 'in1', label='78')
        c79 = Connection(ptpump, 'out1', pt, 'in1', label='79')

        # add connections to network
        self.nw.add_conns(
            c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15,
            c16, c17, c18, c19, c21, c23, c25, c26, c27, c29, c31, c32, c33,
            c34, c35, c36, c37, c39, c40, c41, c42, c44, c45, c46, c47, c49,
            c50, c51, c52, c54, c55, c56, c57, c59, c60, c61, c62, c63, c64,
            c65, c66, c70, c71, c72, c73, c74, c75, c76, c77, c78, c79)

        # power bus
        power = Bus('total output power')
        power.add_comps(
            {'comp': hpt1, 'char': 0.97, 'base': 'component'},
            {'comp': hpt2, 'char': 0.97, 'base': 'component'},
            {'comp': lpt1, 'char': 0.97, 'base': 'component'},
            {'comp': lpt2, 'char': 0.97, 'base': 'component'},
            {'comp': lpt3, 'char': 0.97, 'base': 'component'},
            {'comp': lpt4, 'char': 0.97, 'base': 'component'},
            {'comp': lpt5, 'char': 0.97, 'base': 'component'},
            {'comp': fwp, 'char': 0.95, 'base': 'bus'},
            {'comp': condpump, 'char': 0.95, 'base': 'bus'},
            {'comp': ptpump, 'char': 0.95, 'base': 'bus'},
            {'comp': cwp, 'char': 0.95, 'base': 'bus'},
            {'comp': fan, 'char': 0.95, 'base': 'bus'})

        heat_input_bus = Bus('heat input')
        heat_input_bus.add_comps({'comp': pt, 'base': 'bus'})

        exergy_loss_bus = Bus('exergy loss')
        exergy_loss_bus.add_comps(
            {'comp': air_in, 'base': 'bus'}, {'comp': air_out})

        self.nw.add_busses(power, heat_input_bus, exergy_loss_bus)

        # component parameters
        pt.set_attr(doc=0.95, aoi=0,
                    Tamb=25, A='var', eta_opt=0.73,
                    c_1=0.00496, c_2=0.000691, E=1000,
                    iam_1=1, iam_2=1)

        ptpump.set_attr(eta_s=0.6)

        eco.set_attr()
        eva.set_attr(ttd_l=5)
        sup.set_attr()

        hpt1.set_attr(eta_s=0.8376)
        hpt2.set_attr(eta_s=0.8463)
        lpt1.set_attr(eta_s=0.8623)
        lpt2.set_attr(eta_s=0.917)
        lpt3.set_attr(eta_s=0.9352)
        lpt4.set_attr(eta_s=0.88)
        lpt5.set_attr(eta_s=0.6445)

        cond.set_attr(pr1=1, pr2=0.9, ttd_u=5)
        condpump.set_attr(eta_s=0.7)
        fwp.set_attr(eta_s=0.7)
        cwp.set_attr(eta_s=0.7)
        ct.set_attr(pr1=0.95)
        fan.set_attr(eta_s=0.6)

        lppre1.set_attr(pr1=1, ttd_u=5)
        lppre2.set_attr(pr1=1, ttd_u=5)
        lppre3.set_attr(pr1=1, ttd_u=5)
        hppre1.set_attr(pr1=1, ttd_u=5)
        hppre2.set_attr(pr1=1, ttd_u=5)

        lppre1_sub.set_attr(pr1=1, pr2=1, ttd_l=10)
        lppre2_sub.set_attr(pr1=1, pr2=1, ttd_l=10)
        lppre3_sub.set_attr(pr1=1, pr2=1, ttd_l=10)
        hppre1_sub.set_attr(pr1=1, pr2=1, ttd_l=10)
        hppre2_sub.set_attr(pr1=1, pr2=1, ttd_l=10)

        # connection parameters
        # parabolic trough cycle
        c70.set_attr(fluid={'TVP1': 1, 'water': 0, 'air': 0}, T=390, p=23.304)
        c76.set_attr(m=Ref(c70, 0.1284, 0))
        c73.set_attr(p=22.753)
        c74.set_attr(p=21.167)
        c78.set_attr(p=20.34)
        c79.set_attr(p=41.024)

        # cooling water
        c62.set_attr(
            fluid={'TVP1': 0, 'water': 1, 'air': 0}, T=30, p=self.pamb)
        # cooling tower
        c64.set_attr(
            fluid={'water': 0, 'TVP1': 0, 'air': 1}, p=self.pamb, T=self.Tamb)
        c65.set_attr(p=self.pamb + 0.0005)
        c66.set_attr(p=self.pamb, T=30)
        # power cycle
        c32.set_attr(Td_bp=-2)
        c34.set_attr(x=0.5)
        c1.set_attr(fluid={'water': 1, 'TVP1': 0, 'air': 0}, p=100, T=371)

        # steam generator pressure values
        c31.set_attr(p=103.56)
        c35.set_attr(p=103.42)

        # turbine pressure values
        c3.set_attr(p=33.61, m=38.969)
        c5.set_attr(p=18.58)
        c7.set_attr(p=17.1, T=371)
        c8.set_attr(p=7.98)
        c10.set_attr(p=2.73)
        c12.set_attr(p=0.96)
        c14.set_attr(p=0.29)

        # preheater pressure values
        c19.set_attr(p=14.755, state='l')
        c21.set_attr(p=9.9975, state='l')
        c23.set_attr(p=8.7012, state='l')
        c25.set_attr(state='l')

        c27.set_attr(p=125)
        c29.set_attr(p=112)

        # condensation
        c16.set_attr(p=0.08)

        # feedwater tank
        c26.set_attr(x=0)

        # a stable solution is generated for parts of the network
        self.nw.solve(mode='design')

        self.nw.del_conns(c19, c21, c23, c27, c29, c37, c42, c47, c52, c57)

        c19 = Connection(condpump, 'out1', lppre1_sub, 'in2', label='19')
        c20 = Connection(lppre1_sub, 'out2', lppre1, 'in2', label='20')
        c21 = Connection(lppre1, 'out2', lppre2_sub, 'in2', label='21')
        c22 = Connection(lppre2_sub, 'out2', lppre2, 'in2', label='22')
        c23 = Connection(lppre2, 'out2', lppre3_sub, 'in2', label='23')
        c24 = Connection(lppre3_sub, 'out2', lppre3, 'in2', label='24')

        c27 = Connection(fwp, 'out1', hppre1_sub, 'in2', label='27')
        c28 = Connection(hppre1_sub, 'out2', hppre1, 'in2', label='28')
        c29 = Connection(hppre1, 'out2', hppre2_sub, 'in2', label='29')
        c30 = Connection(hppre2_sub, 'out2', hppre2, 'in2', label='30')

        c37 = Connection(hppre2, 'out1', hppre2_sub, 'in1', label='37')
        c38 = Connection(hppre2_sub, 'out1', v1, 'in1', label='38')
        c42 = Connection(hppre1, 'out1', hppre1_sub, 'in1', label='42')
        c43 = Connection(hppre1_sub, 'out1', v2, 'in1', label='43')

        c47 = Connection(lppre3, 'out1', lppre3_sub, 'in1', label='47')
        c48 = Connection(lppre3_sub, 'out1', v3, 'in1', label='48')
        c52 = Connection(lppre2, 'out1', lppre2_sub, 'in1', label='52')
        c53 = Connection(lppre2_sub, 'out1', v4, 'in1', label='53')
        c57 = Connection(lppre1, 'out1', lppre1_sub, 'in1', label='57')
        c58 = Connection(lppre1_sub, 'out1', v5, 'in1', label='58')

        self.nw.add_conns(
            c19, c20, c21, c22, c23, c24, c27, c28, c29, c30, c37, c38, c42,
            c43, c47, c48, c52, c53, c57, c58)

        # specification of missing parameters
        c19.set_attr(p=14.755)
        c21.set_attr(p=9.9975, state='l')
        c23.set_attr(p=8.7012, state='l')
        c27.set_attr(p=125)
        c29.set_attr(p=112)

        # solve final state
        self.nw.solve(mode='design')

    def test_model(self):
        """Test the thermodynamic model."""
        power_ebsilon = -31.769
        power_tespy = round(
            self.nw.busses['total output power'].P.val / 1e6, 3)
        msg = (
            'The total power calculated (' + str(power_tespy) + ') does not '
            'match the power calculated with the EBSILON model (' +
            str(power_ebsilon) + ').')
        assert power_tespy == power_ebsilon, msg

        T_c79_ebsilon = 296.254
        T_c79_tespy = round(self.nw.get_conn('79').T.val, 3)
        msg = (
            'The temperature at connection 79 calculated (' +
            str(T_c79_tespy) + ') does not match the temperature calculated '
            'with the EBSILON model (' + str(T_c79_ebsilon) + ').')
        assert T_c79_tespy == T_c79_ebsilon, msg

    def test_exergy_analysis(self):
        """Test the exergy analysis results."""
        # carry out exergy analysis
        ean = ExergyAnalysis(
            self.nw, E_P=[self.nw.busses['total output power']],
            E_F=[self.nw.busses['heat input']],
            E_L=[self.nw.busses['exergy loss']])
        ean.analyse(pamb=self.pamb, Tamb=self.Tamb)

        # generate Grassmann diagram
        links, nodes = ean.generate_plotly_sankey_input()

        # check if exergy product value in links is equal to total power
        # output
        position = links['target'].index(nodes.index('E_P'))
        power_links = round(links['value'][position], 0)
        power_bus = round(-self.nw.busses['total output power'].P.val, 0)
        msg = (
            'The exergy product value in the links (' + str(power_links) +
            ') must be equal to the power on the respective bus (' +
            str(power_bus) + ').')
        assert power_links == power_bus, msg
