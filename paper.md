---
title: 'Thermal Engineering Systems in Python'
tags:
  - Thermal Engineering
  - Thermodynamics
  - Power plant simulation
  - Python
authors:
  - name: Francesco Witte
    orcid: 0000-0003-4019-0390
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: Center for Sustainable Energy Systems, Flensburg University of Applied Sciences, Germany
   index: 1
date: someday somemonth 2020
bibliography: paper.bib

---

# Summary

TESPy - Thermal Engineering Systems in Python provides a powerful simulation
toolkit for thermal engineering plants such as power plants, district heating
systems or heat pumps.With the TESPy package it is possible to simulate
stationary operation for process design of thermal energy applications. From
that point it is possible to predict partload behaviour using underlying characteristics for each of the plant's components.

The package inlcudes basic components, such as turbines, pumps, compressors,
heat exchangers, pipes, mixers and splitters as well as some advanced
components. To simulate the operation of a specific power plant an individual
model is created by connecting the respective components to form a topological
network. TESPy automatically generates a set of nonlinear equations based on
the components applied. Using the multi-dimensional Newton-Raphson method the
system of equations is solved to calculate mass flow as well as the fluid
properties pressure and enthalpy along with the fluid composition at every
point of the network.

A CAES power plant is defined by its topology as well as the parametrization of the respective components and connections. The component based power plant simulation software TESPy provides pre-defined components, for example different types of turbomachinery, heat exchangers and combustion chambers. Each component comes with built-in basic equations (e.g. mass flow balance) and optional equations (e.g. specification of isentropic efficiency).
To simulate the operation of a specific power plant an individual model is created by connecting the respective components to form a topological network. TESPy automatically generates a set of nonlinear equations based on the components applied. Using the multi-dimensional Newton-Raphson method the system of equations is solved to calculate mass flow as well as fluid properties such as pressure, enthalpy and fluid composition at every point of the network.
Description of the software

Various technologies with individual topologies allow the integration of geothermal storage sites into an energy supply system. Thus, the power plant model integrated in a coupled simulator must provide usage of generic topologies within the technology scope. TESPy is a component based power plant simulation software that provides a variety of pre-defined components such as different types of turbomachinery, heat exchangers and piping components. To simulate the operation of a specific plant an individual model is created by connecting the respective components to form a topological network. TESPy automatically generates a set of nonlinear equations based on the components applied. Using the multi-dimensional Newton-Raphson method the system of equations is solved to calculate the primary variables mass flow as well as the fluid properties pressure and enthalpy at every point of the network [REFERENCES].
TESPy applies balance equations for every component regarding
    • mass flow,
    • energy and
    • fluid composition.
As the simulation is in steady state mode, the total mass flow into a component must be equal to the total mass flow leaving the component (eq. 1). Additionally, all components except the combustion chamber and merges, do not change the fluid composition, thus equation (2) is applied. The energy balance of all components is derived from the stationary energy balance of open thermodynamic systems with multiple inlets and outlets. Differences in flow velocity as well as height are neglected as these are relatively small compared to change in enthalpy in thermal engineering applications. The values of heat and power transferred, depends on the individual component properties. For example, a pipe does not transfer power, thus only heat may be transferred. In contrast, the turbine is considered adiabatic, thus only transfers power (eq. 3).


The operation of the CAES is calculated with the component based power plant simulation software TESPy. It is designed to simulate stationary operation of thermal power plants and provides a variety of pre-defined components such as different types of turbomachinery, heat exchangers and combustion chambers. Each component comes with built-in basic equations and optional equations.
To simulate the operation of a specific power plant an individual model is created by connecting the respective components to form a topological network. TESPy automatically generates a set of nonlinear equations based on the components applied. Using the multi-dimensional Newton-Raphson method the system of equations is solved to calculate the primary variables mass flow as well as the fluid properties pressure and enthalpy at every point of the network. If multiple fluids are applied within one network, the fluid composition is part of the solution, too.
The following part provides a detailed description of the governing equations for the different components of a CAES in TESPy. Two individual models are created, one for the compression part and one for the expansion part of the power plant. Subsequently, the different CAES setups are presented.
In general, TESPy applies balance equations for every component regarding
    • mass flow,
    • energy and
    • fluid composition.
As the simulation is in steady state mode, the total mass flow into a component must be equal to the total mass flow leaving the component (eq. 1). Additionally, all components except the combustion chamber and merges, do not change the fluid composition, thus equation (2) is applied. The energy balance of all components is derived from the stationary energy balance of open thermodynamic systems with multiple inlets and outlets. Differences in flow velocity as well as height are neglected as these are relatively small compared to change in enthalpy in thermal engineering applications. The values of heat and power transferred, depends on the individual component properties. For example, a pipe does not transfer power, thus only heat may be transferred. In contrast, the turbine is considered adiabatic, thus only transfers power (eq. 3). The energy balance of the combustion chamber is different, as the calorific value of the fuel has to be considered, too. Additionally, it is necessary to compensate for the different zero point definitions for enthalpy in the fluid properties of the fuel, the air and the flue gas. Thus a reference state  is defined.


# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this: ![Example figure.](figure.png)

# Acknowledgements

Acknowledgement to contributors/institution

# References
