# -*- coding: utf-8

"""Module of class TurboCompressor.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/turbomachinery/turbocompressor.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import component_registry
from tespy.components.turbomachinery.compressor import Compressor
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.fluid_properties import isentropic


@component_registry
class TurboCompressor(Compressor):
    r"""
    Class for a turbocompressor.

    **Mandatory Equations**

    - fluid: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - mass flow: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.dp_structure_matrix`
    - :py:meth:`tespy.components.component.Component.pr_structure_matrix`
    - :py:meth:`tespy.components.turbomachinery.base.Turbomachine.energy_balance_func`
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.eta_s_func`
    - :py:meth:`tespy.components.turbomachinery.turbocompressor.TurboCompressor.char_map_eta_s_func`
    - :py:meth:`tespy.components.turbomachinery.turbocompressor.TurboCompressor.char_map_pr_func`

    Inlets/Outlets

    - in1
    - out1

    Optional inlets

    - power

    Image

    .. image:: /api/_images/Compressor.svg
       :alt: flowsheet of the compressor
       :align: center
       :class: only-light

    .. image:: /api/_images/Compressor_darkmode.svg
       :alt: flowsheet of the compressor
       :align: center
       :class: only-dark

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    design_path : str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings : boolean
        Ignore warnings on default characteristics usage for this component.

    printout : boolean
        Include this component in the network's results printout.

    P : float, dict
        Power, :math:`P/\text{W}`

    eta_s : float, dict
        Isentropic efficiency, :math:`\eta_s/1`

    pr : float, dict
        Outlet to inlet pressure ratio, :math:`pr/1`

    dp : float, dict
        Inlet to outlet pressure difference, :math:`dp/\text{p}_\text{unit}`
        Is specified in the Network's pressure unit

    char_map_pr : tespy.tools.characteristics.CharMap, dict
        Characteristic map for pressure ratio vs. nondimensional mass flow.

    char_map_eta_s : tespy.tools.characteristics.CharMap, dict
        Characteristic map for isentropic efficiency vs. nondimensional mass
        flow.

    igva : float, dict, :code:`"var"`
        Inlet guide vane angle, :math:`igva/^\circ`.

    Example
    -------
    Create an air compressor model and calculate the power required for
    compression of 50 l/s of ambient air to 5 bars. Using a generic compressor
    map how does the efficiency change in different operation mode (e.g. 90 %
    of nominal volumetric flow)?

    >>> from tespy.components import Sink, Source, TurboCompressor
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import os
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC", "volumetric_flow": "l/s",
    ...     "enthalpy": "kJ/kg"
    ... })
    >>> si = Sink('sink')
    >>> so = Source('source')
    >>> comp = TurboCompressor('compressor')
    >>> inc = Connection(so, 'out1', comp, 'in1')
    >>> outg = Connection(comp, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    Specify the compressor parameters: nominal efficiency and pressure ratio.
    For offdesign mode the characteristic map is selected instead of the
    isentropic efficiency. For offdesign, the inlet guide vane angle should be
    variable in order to maintain the same pressure ratio at a different
    volumetric flow.

    >>> comp.set_attr(
    ...     pr=5, eta_s=0.8, design=['eta_s'],
    ...     offdesign=['char_map_pr', 'char_map_eta_s']
    ... )
    >>> inc.set_attr(fluid={'air': 1}, p=1, T=20, v=50)
    >>> nw.solve('design')
    >>> nw.save('tmp.json')
    >>> round(comp.P.val, 0)
    12772.0
    >>> round(comp.eta_s.val, 2)
    0.8
    >>> inc.set_attr(v=45)
    >>> comp.set_attr(igva='var')
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(comp.eta_s.val, 2)
    0.77
    >>> os.remove('tmp.json')
    """

    def _preprocess(self, row_idx):
        # skip the FutureWarning of the Compressor class
        return super(Compressor, self)._preprocess(row_idx)

    def get_parameters(self):
        parameters = super().get_parameters()
        parameters.update({
            'igva': dc_cp(min_val=-90, max_val=90, d=1e-4, val=0, quantity="angle"),
            'char_map_eta_s': dc_cm(),
            'char_map_eta_s_group': dc_gcp(
                elements=['char_map_eta_s', 'igva'], num_eq_sets=1,
                func=self.char_map_eta_s_func,
                dependents=self.char_map_dependents
            ),
            'char_map_pr': dc_cm(),
            'char_map_pr_group': dc_gcp(
                elements=['char_map_pr', 'igva'],
                num_eq_sets=1,
                func=self.char_map_pr_func,
                dependents=self.char_map_dependents
            )
        })
        del parameters["eta_s_char"]
        return parameters

    def char_map_pr_func(self):

        r"""
        Calculate pressure ratio from characteristic map.

        Returns
        -------
        residual : float
            Residual value of equations.

        Note
        ----
        - X: speedline index (rotational speed is constant)
        - Y: nondimensional mass flow
        - igva: variable inlet guide vane angle for value manipulation
          according to :cite:`GasTurb2018`.

        .. math::

            X = \sqrt{\frac{T_\mathrm{in,design}}{T_\mathrm{in}}}\\
            Y = \frac{\dot{m}_\mathrm{in} \cdot p_\mathrm{in,design}}
            {\dot{m}_\mathrm{in,design} \cdot p_\mathrm{in} \cdot X}\\
            \vec{Y} = f\left(X,Y\right)\cdot\left(1-\frac{igva}{100}\right)\\
            \vec{Z} = f\left(X,Y\right)\cdot\left(1-\frac{igva}{100}\right)\\
            0 = \frac{p_{out} \cdot p_{in,design}}
            {p_\mathrm{in} \cdot p_\mathrm{out,design}}-
            f\left(Y,\vec{Y},\vec{Z}\right)
        """
        i = self.inl[0]
        o = self.outl[0]

        beta = np.sqrt(i.T.design / i.calc_T())
        y = (i.m.val_SI * i.p.design) / (i.m.design * i.p.val_SI * beta)

        yarr, zarr = self.char_map_pr.char_func.evaluate_x(beta)
        # value manipulation with igva
        yarr *= (1 - self.igva.val_SI / 100)
        zarr *= (1 - self.igva.val_SI / 100)
        pr = self.char_map_pr.char_func.evaluate_y(y, yarr, zarr)

        return (o.p.val_SI / i.p.val_SI) - pr * self.pr.design

    def char_map_eta_s_func(self):
        r"""
        Calculate isentropic efficiency from characteristic map.

        Returns
        -------
        residual : float
            Residual value of equation.

        Note
        ----
        - X: speedline index (rotational speed is constant)
        - Y: nondimensional mass flow
        - igva: variable inlet guide vane angle for value manipulation
          according to :cite:`GasTurb2018`.

        .. math::

            X = \sqrt{\frac{T_\mathrm{in,design}}{T_\mathrm{in}}}\\
            Y = \frac{\dot{m}_\mathrm{in} \cdot p_\mathrm{in,design}}
            {\dot{m}_\mathrm{in,design} \cdot p_\mathrm{in} \cdot X}\\
            \vec{Y} = f\left(X,Y\right)\cdot\left(1-\frac{igva}{100}\right)\\
            \vec{Z}=f\left(X,Y\right)\cdot\left(1-\frac{igva^2}{10000}\right)\\
            0 = \frac{\eta_\mathrm{s}}{\eta_\mathrm{s,design}} -
            f\left(Y,\vec{Y},\vec{Z}\right)
        """
        i = self.inl[0]
        o = self.outl[0]

        x = np.sqrt(i.T.design / i.calc_T())
        y = (i.m.val_SI * i.p.design) / (i.m.design * i.p.val_SI * x)

        yarr, zarr = self.char_map_eta_s.char_func.evaluate_x(x)
        # value manipulation with igva
        yarr *= (1 - self.igva.val_SI / 100)
        zarr *= (1 - self.igva.val_SI ** 2 / 10000)
        eta = self.char_map_eta_s.char_func.evaluate_y(y, yarr, zarr)

        return (
            (
            isentropic(
                i.p.val_SI,
                i.h.val_SI,
                o.p.val_SI,
                i.fluid_data,
                i.mixing_rule,
                T0=i.T.val_SI
            ) - i.h.val_SI)
            / (o.h.val_SI - i.h.val_SI) - eta * self.eta_s.design
        )

    def char_map_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
            self.igva
        ]

    def check_parameter_bounds(self):
        r"""Check parameter value limits."""
        _no_limit_violations = super().check_parameter_bounds()

        for data in [self.char_map_pr, self.char_map_eta_s]:
            if data.is_set:
                x = np.sqrt(self.inl[0].T.design / self.inl[0].T.val_SI)
                y = (
                    (self.inl[0].m.val_SI * self.inl[0].p.design)
                    / (self.inl[0].m.design * self.inl[0].p.val_SI * x)
                )
                yarr = data.char_func.get_domain_errors_x(x, self.label)
                yarr *= (1 - self.igva.val_SI / 100)
                data.char_func.get_domain_errors_y(y, yarr, self.label)

        return _no_limit_violations
