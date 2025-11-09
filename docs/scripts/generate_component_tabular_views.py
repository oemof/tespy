import os
import pandas as pd
from tabulate import tabulate
from tespy.components.component import component_registry
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import GroupedComponentCharacteristics as dc_gcc
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import SimpleDataContainer as dc_simple

def _create_component_instances() -> dict:
    instances = {}
    for cls_name, cls_ref in component_registry.items.items():
        instances[cls_name] = cls_ref(cls_name)

    return instances


def _get_eq_reference(datacontainer):
    if datacontainer.func is None and datacontainer.structure_matrix is None:
        return

    elif datacontainer.func is None:
            basename = datacontainer.structure_matrix.__qualname__
            eq_reference = f"{datacontainer.structure_matrix.__module__}.{basename}"

    else:
        basename = datacontainer.func.__qualname__
        eq_reference = f"{datacontainer.func.__module__}.{basename}"

    return ":py:meth:`" + basename.split(".")[-1] + " <" + eq_reference + ">`"


def collect_component_parameters(instance):
    grouped_parameters = {}
    parameters_with_equation = {}
    characteristic_lines_and_maps = {}
    constraints = {}

    for key, value in instance.get_mandatory_constraints().items():
        eq_reference = _get_eq_reference(value)
        constraints[key] = {
            "eq_reference": eq_reference,
            "description": value.description
        }

    for key, value in instance.parameters.items():
        eq_reference = _get_eq_reference(value)

        if isinstance(value, dc_cp):
            parameters_with_equation[key] = {
                "eq_reference": eq_reference,
                "description": value.description,
                "quantity": value.quantity
            }
        elif isinstance(value, dc_simple):
            parameters_with_equation[key] = {
                "eq_reference": eq_reference,
                "description": value.description,
                "quantity":  None
            }
        elif isinstance(value, dc_cc) or isinstance(value, dc_cm):
            characteristic_lines_and_maps[key] = {
                "eq_reference": eq_reference,
                "description": value.description
            }
        # gcp and gcc are equivalent
        elif isinstance(value, dc_gcp):
            grouped_parameters[key] = {
                "eq_reference": eq_reference,
                "description": value.description,
                "elements": ", ".join([f":code:`{element}`" for element in value.elements])
            }

    return constraints, parameters_with_equation, grouped_parameters, characteristic_lines_and_maps


def collect_component_constraints():
    return


def _indent_block(s: str, spaces: int) -> str:
    pad = " " * spaces
    return s.replace("\n", "\n" + pad)


instances = _create_component_instances()
path = os.path.dirname(os.path.dirname(__file__))
path = os.path.join(path, "modules", "_components_overview.rst")

modules = {}

for cls_name, instance in instances.items():
    # this is required to collect the mandatory constraints
    # it will keep the power balance equation between the component and the
    # (optionally) attached PowerConnection out of the list of equations
    # TODO: Find a better solution
    instance.power_outl = []
    instance.power_inl = []

    cls_ref = instance.__class__
    cls_module = cls_ref.__module__
    parent_module = ".".join(cls_ref.__module__.split(".")[:-1])
    if parent_module not in modules:
        modules[parent_module] = {}

    constraints, parameters, parameter_groups, characteristic_lines = collect_component_parameters(instance)

    modules[parent_module][cls_name] = "\n"
    modules[parent_module][cls_name] += (f".. rubric:: {cls_name}" + "\n" * 2)
    modules[parent_module][cls_name] += (f"Class documentation and example: :py:class:`{cls_name} <{cls_module}.{cls_name}>`")

    outputs = {
        "Table of constraints": constraints,
        "Table of parameters": parameters,
        "Table of parameter groups": parameter_groups,
        "Table of characteristic lines and maps": characteristic_lines,
    }
    headers = {
        "description": "Description",
        "quantity": "Quantity",
        "elements": "Required parameters",
        "eq_reference": "Method"
    }

    for key, value in outputs.items():
        df = pd.DataFrame.from_dict(value).T
        if not df.empty:
            df = df.fillna(":code:`None`")

            modules[parent_module][cls_name] += ("\n" * 2)
            modules[parent_module][cls_name] += (key + "\n" * 2)
            modules[parent_module][cls_name] += (
                tabulate(
                    df[[c for c in headers if c in df.columns]],
                    headers=["Parameter"] + [v for c, v in headers.items() if c in df.columns],
                    tablefmt="rst"
                )
            )

    modules[parent_module][cls_name] += ("\n" * 2)

with open(path, "w", encoding="utf-8") as f:
    f.write(".. tab-set::\n\n")
    for module, data in modules.items():
        f.write(f"    .. tab-item:: {module.split('.')[-1]}\n")
        for cls_info in data.values():
            f.write(_indent_block(cls_info, 8))
            f.write("\n" * 2)
