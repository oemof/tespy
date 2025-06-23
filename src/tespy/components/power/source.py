from tespy.components.component import Component
from tespy.components.component import component_registry


@component_registry
class PowerSource(Component):
    r"""
    A power flow emerges from a PowerSource.

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
    Create a PowerSource and specify a label.

    >>> from tespy.components import PowerSource
    >>> so = PowerSource('a labeled source')
    >>> so.component()
    'power source'
    >>> so.label
    'a labeled source'
    """

    @staticmethod
    def component():
        return 'power source'

    @staticmethod
    def poweroutlets():
        return ["power_out"]

    @staticmethod
    def get_mandatory_constraints():
        return {}
