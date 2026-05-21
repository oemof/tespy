# -*- coding: utf-8

"""On-demand updater for component and connection class docstrings.

Rewrites the middle sections (image block, Ports, Mandatory Equations,
Parameters) of every registered component docstring from live introspection
of :py:meth:`port_schema`, :py:meth:`get_mandatory_constraints`, and
:py:meth:`get_parameters`.  For connection classes only the Parameters section
is regenerated.  The first paragraph(s) and the :code:`Example` section are
preserved verbatim.

Usage
-----
    # Update all registered components and connections in-place
    python docs/scripts/docstring_updater.py

    # Preview without writing files
    python docs/scripts/docstring_updater.py --dry-run

    # Update specific classes only
    python docs/scripts/docstring_updater.py Turbine Compressor --dry-run

    # Update connections only
    python docs/scripts/docstring_updater.py --connections-only

    # Update components only
    python docs/scripts/docstring_updater.py --components-only
"""

import ast
import inspect
import re
import subprocess
import textwrap

from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import FluidComposition as dc_flu
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.data_containers import GroupedComponentCharacteristics as dc_gcc
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import ReferencedFluidProperties as dc_ref
from tespy.tools.data_containers import SimpleDataContainer as dc_simple

# ---------------------------------------------------------------------------
# Base parameters shared by every component
# ---------------------------------------------------------------------------

_BASE_PARAMETERS = [
    ("label", "str", "The label of the component."),
    ("design", "list", "List containing design parameters (stated as String)."),
    ("offdesign", "list", "List containing offdesign parameters (stated as String)."),
    ("design_path", "str", "Path to the components design case."),
    ("local_offdesign", "bool",
     "Treat this component in offdesign mode in a design calculation."),
    ("local_design", "bool",
     "Treat this component in design mode in an offdesign calculation."),
    ("char_warnings", "bool",
     "Ignore warnings on default characteristics usage for this component."),
    ("printout", "bool", "Include this component in the network's results printout."),
]

# ---------------------------------------------------------------------------
# Base parameters shared by every connection
# ---------------------------------------------------------------------------

_CONNECTION_BASE_PARAMETERS = [
    ("label", "str", "The label of the connection."),
    ("design", "list", "List containing design parameters (stated as String)."),
    ("offdesign", "list", "List containing offdesign parameters (stated as String)."),
    ("design_path", "str", "Path to the individual design case for this connection."),
    ("local_offdesign", "bool",
     "Treat this connection in offdesign mode in a design calculation."),
    ("local_design", "bool",
     "Treat this connection in design mode in an offdesign calculation."),
    ("printout", "bool", "Include this connection in the network's results printout."),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _eq_reference(dc):
    """Return a :code:`:py:meth:` cross-reference string for the equation, or :code:`None`."""
    func = getattr(dc, "func", None)
    mat = getattr(dc, "structure_matrix", None)
    target = func if func is not None else mat
    if target is None:
        return None
    qualname = target.__qualname__
    module = target.__module__
    short = qualname.split(".")[-1]
    return f":py:meth:`{short} <{module}.{qualname}>`"


def _collect_constraints(instance):
    """Return :code:`(unconditional, conditional)` mandatory-constraint dicts.

    Conditional constraints are those that only activate when a power or heat
    connector is attached to the component.
    """
    for attr in ("power_outl", "power_inl", "heat_outl", "heat_inl"):
        setattr(instance, attr, [])

    base = instance.get_mandatory_constraints()

    sentinel = [object()]
    for attr in ("power_outl", "power_inl", "heat_outl", "heat_inl"):
        setattr(instance, attr, sentinel)

    with_connectors = instance.get_mandatory_constraints()

    for attr in ("power_outl", "power_inl", "heat_outl", "heat_inl"):
        setattr(instance, attr, [])

    unconditional = {}
    conditional = {}
    for key, value in with_connectors.items():
        entry = {
            "ref": _eq_reference(value),
            "description": getattr(value, "description", None) or key,
        }
        (unconditional if key in base else conditional)[key] = entry

    return unconditional, conditional


def _dc_type(dc):
    """Return the NumPy-style type annotation string for a data container."""
    import math
    if isinstance(dc, dc_cp):
        base = "float, dict"
        return base + ', :code:`"var"`' if getattr(dc, "_allows_var", False) else base
    if isinstance(dc, dc_cc):
        return "tespy.tools.characteristics.CharLine, dict"
    if isinstance(dc, dc_cm):
        return "tespy.tools.characteristics.CharMap, dict"
    if isinstance(dc, (dc_gcp, dc_gcc)):
        return type(dc).__name__
    if isinstance(dc, dc_ref):
        return "Ref"
    if isinstance(dc, dc_prop):
        return "float, Ref"
    if isinstance(dc, dc_flu):
        return "dict"
    if isinstance(dc, dc_simple):
        dtype = getattr(dc, "dtype", None)
        if dtype:
            return dtype
        val = getattr(dc, "val", None)
        if isinstance(val, bool):
            return "bool"
        if isinstance(val, int):
            return "int"
        if isinstance(val, float) and not math.isnan(val):
            return "int, float"
        if isinstance(val, str):
            return "str"
        return "any"
    return "any"


# ---------------------------------------------------------------------------
# Section formatters
# ---------------------------------------------------------------------------

_PORT_LABELS = [
    ("inlets",       "Fluid inlets"),
    ("outlets",      "Fluid outlets"),
    ("powerinlets",  "Power inlets"),
    ("poweroutlets", "Power outlets"),
    ("heatinlets",   "Heat inlets"),
    ("heatoutlets",  "Heat outlets"),
]


def _ports_section(cls):
    schema = cls.port_schema()
    lines = []
    for key, label in _PORT_LABELS:
        entry = schema.get(key, {})
        kind = entry.get("type")
        if kind == "fixed":
            ports = entry.get("ports", [])
            if ports:
                lines.append(f"{label}: {', '.join(ports)}")
        elif kind == "variable":
            param = entry["parameter"]
            pattern = entry["pattern"]
            lines.append(
                f"{label}: {pattern.replace('{n}', '1')}, "
                f"{pattern.replace('{n}', '2')}, ... "
                f"(variable, count set by :code:`{param}`)"
            )
    if not lines:
        return ""
    return "Ports\n-----\n\n" + "\n\n".join(lines)


def _mandatory_section(instance):
    try:
        unconditional, conditional = _collect_constraints(instance)
    except Exception:
        return ""

    if not unconditional and not conditional:
        return "Mandatory Equations\n-------------------\n\nNone"

    def _bullet(entry):
        ref = entry["ref"]
        desc = entry["description"].rstrip(".")
        return "- " + desc + (f": {ref}" if ref else "")

    parts = ["Mandatory Equations\n-------------------\n"]
    parts.extend(_bullet(e) for e in unconditional.values())

    if conditional:
        parts.append("\nWhen a power or heat connector is attached:\n")
        parts.extend(_bullet(e) for e in conditional.values())

    return "\n".join(parts)


def _parameters_section(instance, base_params=None, param_filter=None):
    if base_params is None:
        base_params = _BASE_PARAMETERS
    entries = []

    for name, ptype, desc in base_params:
        entries.append((name, f"\n{name} : {ptype}\n    {desc}"))

    try:
        params = instance.get_parameters()
    except Exception:
        params = {}

    for name, dc in params.items():
        if param_filter is not None and not param_filter(name, dc):
            continue
        ptype = _dc_type(dc)

        raw_desc = getattr(dc, "description", None) or ""
        desc = (raw_desc[0].upper() + raw_desc[1:]).rstrip(".") if raw_desc else ""

        meta_parts = []
        if isinstance(dc, dc_cp):
            q = getattr(dc, "quantity", None)
            if q:
                meta_parts.append(f"Quantity: :code:`{q}`.")
            if getattr(dc, "_allows_var", False):
                meta_parts.append(
                    'Can be set as a system variable by passing :code:`"var"` as its value.'
                )
        if isinstance(dc, (dc_gcp, dc_gcc)):
            elements = getattr(dc, "elements", None) or []
            if elements:
                meta_parts.append(
                    "Elements: " + ", ".join(f":code:`{e}`" for e in elements) + "."
                )

        ref = _eq_reference(dc)

        desc_body = (desc + "." if desc else "") + (
            " " + " ".join(meta_parts) if meta_parts else ""
        )
        body_line = textwrap.fill(
            desc_body.strip(), width=76, initial_indent="    ",
            subsequent_indent="    "
        )
        if ref:
            eq_line = f"    Equation: {ref}."
            entries.append((name, f"\n{name} : {ptype}\n{body_line}\n{eq_line}"))
        else:
            entries.append((name, f"\n{name} : {ptype}\n{body_line}"))

    entries.sort(key=lambda x: x[0].lower())

    return "\n".join(["Parameters\n----------"] + [e for _, e in entries])


def _image_block(cls):
    import os
    _images_dir = os.path.join(os.path.dirname(__file__), "..", "api", "_images", "components")
    available = {
        os.path.splitext(f)[0]
        for f in os.listdir(_images_dir)
        if f.endswith(".svg") and not f.endswith("_darkmode.svg")
    }

    image_name = next(
        (c.__name__ for c in cls.__mro__ if c.__name__ in available),
        None,
    )
    if not image_name:
        return ""

    low = cls.__name__.lower()
    return (
        f".. image:: /api/_images/components/{image_name}.svg\n"
        f"   :alt: flowsheet of the {low}\n"
        f"   :align: center\n"
        f"   :class: only-light\n\n"
        f".. image:: /api/_images/components/{image_name}_darkmode.svg\n"
        f"   :alt: flowsheet of the {low}\n"
        f"   :align: center\n"
        f"   :class: only-dark"
    )


# ---------------------------------------------------------------------------
# Existing-docstring splitting
# ---------------------------------------------------------------------------

# Section markers that signal the start of auto-generated content
_SECTION_START = re.compile(
    r"(?:^|\n)(?=\*\*(?:Mandatory|Optional)"
    r"|Ports\n[-]"
    r"|Mandatory\s+Equations"
    r"|Optional\s+Equations"
    r"|Inlets/Outlets"
    r"|Image\n"
    r"|\.\.\ image::"
    r"|Parameters\n[-])",
    re.MULTILINE,
)
_EXAMPLE_START = re.compile(r"(?:^|\n\n)Examples?\n[-]+", re.MULTILINE)
_PRESERVE_SECTION_START = re.compile(
    r"(?:^|\n\n)(?:"
    r"(?:Note|Notes|Warning|Warnings|References|See Also)\n[-]+"
    r"|\.\.\ (?:note|tip|warning|attention|important|caution)::"
    r")",
    re.MULTILINE,
)


def _split_docstring(cls):
    """Return :code:`(intro, suffix, example)` from the current class docstring.

    *suffix* contains any Note / Warning / References sections that follow
    the auto-generated block; they are preserved verbatim after Parameters.
    """
    doc = inspect.getdoc(cls) or ""

    sm = _SECTION_START.search(doc)
    em = _EXAMPLE_START.search(doc)

    # intro ends at the first auto-generated section or the Example section
    intro_end = len(doc)
    if sm:
        intro_end = sm.start() + (1 if doc[sm.start()] == "\n" else 0)
    if em:
        ex_pos = em.start() + (1 if doc[em.start()] == "\n" else 0)
        intro_end = min(intro_end, ex_pos)

    intro = doc[:intro_end].strip()
    example = doc[em.start():].strip() if em else ""
    if example.startswith("\n"):
        example = example[1:]

    auto_end = em.start() if em else len(doc)
    auto_block = doc[intro_end:auto_end]
    pm = _PRESERVE_SECTION_START.search(auto_block)
    suffix = auto_block[pm.start():].lstrip("\n").strip() if pm else ""

    return intro, suffix, example


# ---------------------------------------------------------------------------
# Docstring assembly
# ---------------------------------------------------------------------------

def generate_docstring(cls):
    """Return the complete new docstring body (without surrounding quotes)."""
    from tespy.components.component import Component
    from tespy.tools.schema import _instantiate_component

    intro, suffix, example = _split_docstring(cls)

    instance = _instantiate_component(cls)
    if instance is None:
        instance = cls.__new__(cls)
        Component.__init__(instance, "_docprobe_")

    parts = []
    if intro:
        parts.append(intro)
    image_block = _image_block(cls)
    if image_block:
        parts.append(image_block)
    ports = _ports_section(cls)
    if ports:
        parts.append(ports)
    mandatory = _mandatory_section(instance)
    if mandatory:
        parts.append(mandatory)
    parts.append(_parameters_section(instance))
    if suffix:
        if re.match(r"^\.\.", suffix):
            parts.append("Notes\n-----\n\n" + suffix)
        else:
            parts.append(suffix)
    if example:
        parts.append(example)

    return "\n\n".join(parts)


# ---------------------------------------------------------------------------
# Source-file patching
# ---------------------------------------------------------------------------

def _char_pos(source, line_1based, col):
    """Convert 1-based line + 0-based col to a character offset in *source*."""
    pos = 0
    for _ in range(line_1based - 1):
        pos = source.index("\n", pos) + 1
    return pos + col


def _patch_file(source, cls_name, new_body):
    """Return *source* with the docstring of *cls_name* replaced by *new_body*."""
    tree = ast.parse(source)

    for node in ast.walk(tree):
        if not (isinstance(node, ast.ClassDef) and node.name == cls_name):
            continue
        if not node.body:
            break
        first = node.body[0]
        if not (
            isinstance(first, ast.Expr)
            and isinstance(first.value, ast.Constant)
            and isinstance(first.value.value, str)
        ):
            break

        dv = first.value
        start = _char_pos(source, dv.lineno, dv.col_offset)
        end = _char_pos(source, dv.end_lineno, dv.end_col_offset)

        indent = " " * dv.col_offset
        raw_prefix = "r" if source[start] == "r" else ""

        # Indent every line of the body to match the class indentation
        indented_lines = []
        for line in new_body.split("\n"):
            indented_lines.append(indent + line if line.strip() else "")
        body = "\n".join(indented_lines)

        new_doc = f'{raw_prefix}"""\n{body}\n{indent}"""'
        return source[:start] + new_doc + source[end:]

    return source


# ---------------------------------------------------------------------------
# Git helpers
# ---------------------------------------------------------------------------

def _has_uncommitted_changes(filepath):
    """Return True if *filepath* has uncommitted changes according to git."""
    try:
        result = subprocess.run(
            ["git", "status", "--porcelain", filepath],
            capture_output=True, text=True, check=True,
        )
        return bool(result.stdout.strip())
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


# ---------------------------------------------------------------------------
# Connection docstring generation
# ---------------------------------------------------------------------------

def generate_connection_docstring(cls):
    """Return the complete new docstring body for a connection class."""
    intro, suffix, example = _split_docstring(cls)

    instance = cls.__new__(cls)

    parts = []
    if intro:
        parts.append(intro)
    parts.append(_parameters_section(
        instance,
        base_params=_CONNECTION_BASE_PARAMETERS,
        param_filter=lambda name, _: not name.endswith("_ref"),
    ))
    if suffix:
        if re.match(r"^\.\.", suffix):
            parts.append("Notes\n-----\n\n" + suffix)
        else:
            parts.append(suffix)
    if example:
        parts.append(example)

    return "\n\n".join(parts)


def update_connection_docstrings(classes=None, dry_run=False, force=False):
    """Update docstrings for *classes* (defaults to all registered connections).

    Parameters
    ----------
    classes : list of str, optional
        Class names to process.  Defaults to every entry in the connection
        registry.
    dry_run : bool
        When :code:`True`, print the generated docstrings to stdout instead of
        writing files.
    force : bool
        When :code:`True`, write files even if they have uncommitted changes.
    """
    from tespy.connections.connection import connection_registry

    registry = {
        k: v for k, v in connection_registry.items.items()
        if classes is None or k in classes
    }

    pending = {}

    for cls_name, cls in registry.items():
        try:
            new_body = generate_connection_docstring(cls)
        except Exception as exc:
            print(f"[SKIP] {cls_name}: {exc}")
            continue

        if dry_run:
            sep = "=" * 72
            print(f"\n{sep}\n=== {cls_name}\n{sep}")
            print(new_body)
            continue

        filepath = inspect.getfile(cls)
        if filepath not in pending:
            with open(filepath, encoding="utf-8") as fh:
                pending[filepath] = fh.read()

        pending[filepath] = _patch_file(pending[filepath], cls_name, new_body)
        print(f"[OK] {cls_name}")

    for filepath, content in pending.items():
        if not force and _has_uncommitted_changes(filepath):
            print(f"[SKIP] {filepath}: has uncommitted changes (pass --force to override)")
            continue
        with open(filepath, "w", encoding="utf-8") as fh:
            fh.write(content)
        print(f"[WRITE] {filepath}")


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def update_component_docstrings(classes=None, dry_run=False, force=False):
    """Update docstrings for *classes* (defaults to all registered components).

    Parameters
    ----------
    classes : list of str, optional
        Class names to process.  Defaults to every entry in the component
        registry except the base :code:`Component` class.
    dry_run : bool
        When :code:`True`, print the generated docstrings to stdout instead of
        writing files.
    force : bool
        When :code:`True`, write files even if they have uncommitted changes.
    """
    from tespy.components.component import component_registry

    registry = {
        k: v for k, v in component_registry.items.items()
        if k != "Component" and (classes is None or k in classes)
    }

    pending = {}  # filepath → current source text

    for cls_name, cls in registry.items():
        try:
            new_body = generate_docstring(cls)
        except Exception as exc:
            print(f"[SKIP] {cls_name}: {exc}")
            continue

        if dry_run:
            sep = "=" * 72
            print(f"\n{sep}\n=== {cls_name}\n{sep}")
            print(new_body)
            continue

        filepath = inspect.getfile(cls)
        if filepath not in pending:
            with open(filepath, encoding="utf-8") as fh:
                pending[filepath] = fh.read()

        pending[filepath] = _patch_file(pending[filepath], cls_name, new_body)
        print(f"[OK] {cls_name}")

    for filepath, content in pending.items():
        if not force and _has_uncommitted_changes(filepath):
            print(f"[SKIP] {filepath}: has uncommitted changes (pass --force to override)")
            continue
        with open(filepath, "w", encoding="utf-8") as fh:
            fh.write(content)
        print(f"[WRITE] {filepath}")


if __name__ == "__main__":
    import sys

    argv = sys.argv[1:]
    dry_run = "--dry-run" in argv
    force = "--force" in argv
    components_only = "--components-only" in argv
    connections_only = "--connections-only" in argv
    names = [a for a in argv if not a.startswith("--")] or None

    if not connections_only:
        update_component_docstrings(classes=names, dry_run=dry_run, force=force)
    if not components_only:
        update_connection_docstrings(classes=names, dry_run=dry_run, force=force)
