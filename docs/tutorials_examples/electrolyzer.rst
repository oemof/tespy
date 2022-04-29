Water Electrolyzer
---------------

This example shows how you can use the water electrolyzer component in TESPy.

.. figure:: _images/WaterElectrolyzer.svg
    :align: center

    Figure: Topology of the water electrolyzer.

The water electrolyzer produces hydrogen and oxygen from water and power.

Example
#######

Create a water electrolyzer and compress the hydrogen, e.g. for a hydrogen storage.

.. code-block:: python

    >>> from tespy.components import (Sink, Source, Compressor, WaterElectrolyzer)
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import shutil
    >>> import numpy as np
    >>> fluid_list = ['O2', 'H2O', 'H2']
    >>> nw = Network(fluids=fluid_list, T_unit='C', p_unit='bar',
    ... v_unit='l / s', iterinfo=False)
    >>> fw = Source('feed water')
    >>> oxy = Sink('oxygen sink')
    >>> hydro = Sink('hydrogen sink')
    >>> cw_cold = Source('cooling water source')
    >>> cw_hot = Sink('cooling water sink')
    >>> comp = Compressor('compressor', eta_s=0.9)
    >>> el = WaterElectrolyzer('electrolyzer')
    >>> el.component()
    'water electrolyzer'

The electrolyzer should produce 100 l/s of hydrogen at an operating pressure of 10 bars and an outlet temperature of 50 °C. The fluid composition needs to be specified for the cooling liquid only. The storage pressure is 25 bars. The electrolysis efficiency is at 80 % and the compressor isentropic efficiency at 85 %. After designing the plant the offdesign electrolysis efficiency is predicted by the characteristic line. The default characteristic line can be found here: :py:mod:`tespy.data`.

.. code-block:: python

    >>> fw_el = Connection(fw, 'out1', el, 'in2')
    >>> el_o = Connection(el, 'out2', oxy, 'in1')
    >>> el_cmp = Connection(el, 'out3', comp, 'in1')
    >>> cmp_h = Connection(comp, 'out1', hydro, 'in1')
    >>> cw_el = Connection(cw_cold, 'out1', el, 'in1')
    >>> el_cw = Connection(el, 'out1', cw_hot, 'in1')
    >>> nw.add_conns(fw_el, el_o, el_cmp, cmp_h, cw_el, el_cw)
    >>> fw_el.set_attr(p=10, T=15)
    >>> cw_el.set_attr(p=5, T=15, fluid={'H2O': 1, 'H2': 0, 'O2': 0})
    >>> el_cw.set_attr(T=45)
    >>> cmp_h.set_attr(p=25)
    >>> el_cmp.set_attr(v=100, T=50)
    >>> el.set_attr(eta=0.8, pr=0.99, design=['eta', 'pr'],
    ... offdesign=['eta_char', 'zeta'])
    >>> comp.set_attr(eta_s=0.85)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(el.e0 / el.P.val * el_cmp.m.val_SI, 1)
    0.8
    >>> P_design = el.P.val / 1e6
    >>> round(P_design, 1)
    13.2
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(el.eta.val, 1)
    0.8
    >>> el_cmp.set_attr(v=np.nan)
    >>> el.set_attr(P=P_design * 0.66)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(el.eta.val, 2)
    0.88
    >>> shutil.rmtree('./tmp', ignore_errors=True)
