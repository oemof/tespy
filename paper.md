---
title: 'TESPy: Thermal Engineering Systems in Python'

tags:
  - Thermal Engineering
  - Thermodynamics
  - Power Plant Simulation
  - Python

authors:
  - name: Francesco Witte
    orcid: 0000-0003-4019-0390
    affiliation: "1, 2"
  - name: Ilja Tuschy
    affiliation: "1, 2"

affiliations:
 - name: Center for Sustainable Energy Systems, Flensburg
   index: 1
 - name: Flensburg University of Applied Sciences
   index: 2

date: 18. February 2020

bibliography: paper.bib

---

# Summary

TESPy (Thermal Engineering Systems in Python) provides a powerful
simulation toolkit for thermal process engineering, for instance power plants,
district heating systems or heat pumps. With the TESPy package it is possible to
design plants and simulate stationary operation. From that point, part-load
behavior can be predicted using underlying characteristics for each of the
plant's components. The component based structure combined with the solution
method offer very high flexibility regarding the plant's topology and its
parametrisation.

# Motivation

Thermal process simulation is a fundamental discipline in energy engineering.
Consequently, there are several well known commercial software solutions
available in this field, for example EBSILON Professional or Aspen Plus, mainly
used in the industrial environment. In order to open this kind of software to a
wide (scientific) audience, an open source solution is very promising. Its
prevalence in the scientific field and the availability of interfaces to other
programming languages plead for a python implementation.

# Method

The package is built of basic components, for example different types of
turbomachinery, heat exchangers, pipes, valves and mixers or splitters. To
simulate a specific plant an individual model is created connecting the
components to form a topological network. Steady state operation of the plant is
determined by the fluid's state on every connection (and therefore its change
within components) between the plant's individual components. Thus, based on the
components and parameters applied, TESPy generates a set of nonlinear equations
representing the plant's topology and its parametrisation. By formulating the
equations implicitly parameters and results are generally interchangeable
offering high flexibility in the user specifications. The system of equations is
numerically solved with an inbuilt solver applying the multi-dimensional
Newton-Raphson method to determine the mass flow, pressure, enthalpy and fluid
composition of every connection. This way, it is possible to solve for both
thermal as well as hydraulic state of the plant.

To achieve this, TESPy implements balance equations (based on standard
literature, for example [@Baehr; @Epple]) for every component regarding

• mass flow and fluid composition as well as

• energy (covering thermal and hydraulic properties).

In steady state, the total mass flow into a component must be equal to the total
mass flow leaving the component (eq. 1). Additionally, the mass balance of a
specific fluid fl is applied (eq. 2). The energy balance of all components is
derived from the stationary energy balance of open systems with multiple inlets
and outlets. Differences in flow velocity as well as height are neglected as
these are relatively small compared to change in enthalpy in thermal engineering
applications. The values of heat and power transferred, depend on the individual
component properties. For example, a pipe does not transfer power, thus only
heat may be transferred. In contrast, turbomachinery is considered adiabatic,
thus only transfers power (eq. 3). If chemical reactions take place, the
corresponding chemical mass balance is taken into account instead of equation 2.
On top of that, the energy balance is different, as the reaction enthalpy has to
be considered, too. Furthermore, it is necessary to compensate for the different
zero point definitions of enthalpy in the fluid properties of the reaction
components by defining a reference state. E. g. for a combustion chamber
equation 4 is implemented.

$$0=\underset{i}{\sum}\dot{m}_{\mathrm{in,}i}-\underset{o}{\sum}\dot{m}_{\mathrm{out,}o}$$

$$0=\underset{i}{\sum}\dot{m}_{\mathrm{in,}i}\cdot x_{fl\mathrm{,in,}i}-
\underset{o}{\sum}\dot{m}_{\mathrm{out,}o}\cdot x_{fl\mathrm{,out,}o}
\ \forall fl\in\mathrm{network\ fluids}$$

$$0=\underset{o}{\sum}\dot{m}_{\mathrm{out,}o}\cdot h_{\mathrm{out,}o}-
\underset{i}{\sum}\dot{m}_{\mathrm{in,}i}\cdot h_{\mathrm{in,}i}-P-\dot{Q}$$

$$0=
\dot{m}_{\mathrm{out}}\cdot\left(h_{\mathrm{out}}-h_{\mathrm{out,ref}}\right)-
\underset{i}{\sum}\dot{m}_{\mathrm{in,}i}\cdot\left(h_{\mathrm{in}}-
h_{\mathrm{in,ref}}\right)_{i}-
\dot{m}_{\mathrm{in}}\cdot\underset{j}{\sum}LHV_{j}\cdot x_{j}$$

After designing a specific plant, part-load performance can be determined. For
this, design specific component parameters are calculated in the design case,
for example the area independent heat transfer coefficient $kA$ of heat
exchangers. The heat transferred at a different operation point may then be
calculated using the design value of $kA$ applying equation (5).

$$0=\dot{Q}-kA\cdot\Delta\vartheta_{\mathrm{log}}$$

In general, these parameters can be adjusted using lookup table functions to
match the model behavior to measured data. This is especially useful if
components with well known characteristics should be implemented in a different
plant or at different operating conditions. Due to the modular structure of
TESPy, new equations or characteristic functions to further improve the
representation of an actual plant can easily be added.

# Example Implementations

The core strength of TESPy lies in the generic and component based architecture
allowing to simulate technologically and topologically different thermal
engineering applications.

For example, as the increasing share of renewable energy sources to mitigate
climate change will result in a significant storage demand, underground gas
storage is considered a large scale energy storage option [@IEA; @Pfeiffer]. Due
to the feedback regarding the physical parameters of the fluid exchanged between
the geological storage and the above-ground plant, an integrated simulation of
the storage and the power plant is necessary for a detailed assessment of such
storage options [@CAES]. Another important task in energy system transition is
renewable heating: Heat pumps using subsurface heat to provide heating on
household or even district level show analogous feedback reactions with the
underground. As their electricity consumption highly depends on the heat source
temperature level, simulator coupling provides valuable assessment possibilities
in this field, too. Additionally, TESPy has been coupled with OpenGeoSys
[@ogs] for pipeline network simulation of borehole thermal energy storage arrays
[@BTES].

# Acknowledgements

This work is supported by University of Applied Sciences and the Center for
Sustainable Energy Systems in Flensburg. It is part of the open energy modeling
framework (oemof) [@oemof]. Many thanks to all
[contributers](https://github.com/oemof/tespy/graphs/contributors).

Key parts of TESPy require the following scientific software packages: CoolProp
[@CoolProp], NumPy [@NumPy], pandas [@pandas]. Other packages implemented are
tabulate and SciPy.

# References
