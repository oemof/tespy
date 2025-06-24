from tespy.components.component import Component
from tespy.components.component import component_registry


@component_registry
class PowerSink(Component):
    r"""
    A power flow drains in a PowerSink.

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

    Example
    -------
    Create a PowerSink and specify a label.

    >>> from tespy.components import PowerSink
    >>> si = PowerSink('a labeled sink')
    >>> si.label
    'a labeled sink'
    """

    @staticmethod
    def powerinlets():
        return ["power"]

    @staticmethod
    def get_mandatory_constraints():
        return {}
