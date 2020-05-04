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
used in the industrial environment. A similar approach was carried out by
[@ThermoPower] with the open source library "ThermoPower" written in the
Modelica language. The software provides dynamic modeling and focuses on steam
power plants, but has not been updated since 2014 [@ThermoPowerOnline].

Due to the fact that ThermoPower is not updated anymore and limited on
one field of thermal engineering a new solution is proposed. In order to open
the software to a wide (scientific) audience, an open source solution
implemented in a widespread programming language like Python, is very promising.

# Package Architecture

TESPy is built up in modular structure with three main modules:

* networks: The networks module represents the container for every simulation.
  The solution process, pre- and postprocessing are handled by the network
  class.
* components: All components are part of the components module. The package
  provides basic components, for example different types of turbomachinery, heat
  exchangers, pipes, valves and mixers or splitters. On top of that, advanced
  components like a drum, an internal combustion engine or a combustion chamber
  are available. The individual properties of each component are defined in the respective class.
* connections: Connections are used to connect components and hold the fluid
  property information. Thus, they represent the plant's topology and its state.

# Method

To simulate a specific plant an individual model is created connecting the
components to form a topological network. Steady state operation of the plant is
determined by the fluid's state on every connection (and therefore its change
within components) between the plant's individual components. Thus, based on the
components and parameters applied, TESPy generates a set of nonlinear equations
representing the plant's topology and its parametrisation. By formulating the
equations implicitly parameters and results are generally interchangeable
offering high flexibility in the user specifications. The system of equations is
numerically solved with an inbuilt solver applying the multi-dimensional
Newton-Raphson method to determine the mass flow $\dot{m}$, pressure $p$,
enthalpy $h$ and fluid composition (defined by mass fraction $x$ of each fluid)
of every connection. This way, it is possible to solve for both thermal as well
as hydraulic state of the plant.

To achieve this, TESPy implements balance equations (based on standard
literature, for example [@Baehr; @Epple]) for every component regarding

* mass flow and fluid composition as well as
* energy (covering thermal and hydraulic properties).

In steady state, the total mass flow into a component must be equal to the total
mass flow leaving the component (eq. \ref{eq:1}). Additionally, the mass balance
of a specific fluid $fl$ is applied (eq. \ref{eq:2}). The energy balance of all
components is derived from the stationary energy balance of open systems with
multiple inlets and outlets (eq. \ref{eq:3}).

In thermal engineering applications the change in kinetic and potential energy
due to differences in flow velocity $c$ and height $z$ are usually neglected as
these are relatively small compared to change in enthalpy, subsequently the
equation can be simplified (eq. \ref{eq:4}).

\begin{align}
\begin{split}
0 &=\underset{i}{\sum}\dot{m}_{\mathrm{in,}i}-\underset{o}{\sum}
\dot{m}_{\mathrm{out,}o} \label{eq:1}
\end{split}
\\[2ex]
\begin{split}
0 &=\underset{i}{\sum}\dot{m}_{\mathrm{in,}i}\cdot x_{fl\mathrm{,in,}i}-
\underset{o}{\sum}\dot{m}_{\mathrm{out,}o}\cdot x_{fl\mathrm{,out,}o}
\ \forall fl\in\mathrm{network\ fluids} \label{eq:2}
\end{split}
\\[2ex]
\begin{split}
0 & =\underset{o}{\sum}\dot{m}_{\mathrm{out,}o}\cdot \left(
h_{\mathrm{out,}o} + \mathrm{g} \cdot z_{\mathrm{out,}o} +
\frac{c_{\mathrm{out,}o}^2}{2}\right)\\ &-
\underset{i}{\sum}\dot{m}_{\mathrm{in,}i}\cdot \left(
h_{\mathrm{in,}i} + \mathrm{g} \cdot z_{\mathrm{in,}i} +
\frac{c_{\mathrm{in,}i}^2}{2}\right)-P-\dot{Q} \label{eq:3}
\end{split}
\\[2ex]
\begin{split}
0 &=\underset{o}{\sum}\dot{m}_{\mathrm{out,}o}\cdot h_{\mathrm{out,}o}-
\underset{i}{\sum}\dot{m}_{\mathrm{in,}i}\cdot h_{\mathrm{in,}i}-P-\dot{Q}
\label{eq:4}
\end{split}
\end{align}

The values of heat $\dot{Q}$ and power $P$ transferred, depend on
the individual component properties. For example,

* a pipe does not transfer power, thus only heat may be transferred ($P=0$).
* turbomachinery is considered adiabatic, thus it only transfers power
  ($\dot{Q}=0$).

If chemical reactions take place, the corresponding chemical mass balance is
taken into account instead of equation \ref{eq:2}. On top of that, the energy
balance is different, as the reaction enthalpy has to be considered.
Furthermore, it is necessary to compensate for the different zero point
definitions of enthalpy in the fluid properties of the reaction components by
defining a reference state $\mathrm{ref}$. E.g. for an adiabatic combustion
chamber equation \ref{eq:5} is implemented using the fuel's lower heating value
$LHV$ as reaction enthalpy.

\begin{equation}
0=
\dot{m}_{\mathrm{out}}\cdot\left(h_{\mathrm{out}}-h_{\mathrm{out,ref}}\right)-
\underset{i}{\sum}\dot{m}_{\mathrm{in,}i}\cdot\left(h_{\mathrm{in}}-
h_{\mathrm{in,ref}}\right)_{i}-
\dot{m}_{\mathrm{in}}\cdot\underset{j}{\sum}LHV_{j}\cdot x_{j} \label{eq:5}
\end{equation}

With respect to hydraulic state, the pressure drop
$p_\mathrm{out}-p_\mathrm{in}$ in a pipe can be calculated
by the Darcy-Weisbach equation \ref{eq:6} from its dimensions (length $L$,
diameter $D$) and the Darcy friction factor $\lambda$ calculated from the
pipe's roughness $k_\mathrm{s}$, its diameter and the Reynolds number $re$.
Another much more simple approach is to specify the pressure ratio $pr$
\ref{eq:7}.

\begin{align}
0=p_\mathrm{out}-p_\mathrm{in} +
\frac{\rho \cdot c^2 \cdot \lambda\left(re,k_\mathrm{s},D\right) \cdot L}
{2 \cdot D} \label{eq:6}\\
0=pr - \frac{p_\mathrm{out}}{p_\mathrm{in}} \label{eq:7}
\end{align}

After designing a specific plant, part-load performance can be determined. For
this, design specific component parameters are calculated in the design case,
for example: the area independent heat transfer coefficient $kA$ of heat
exchangers. The heat transfer at a different operation point is calculated from
the $kA$ value and the logarithmic temperature difference
$\Delta\vartheta_\mathrm{log}$ in equation \ref{eq:8}.

 \begin{equation}
 0=\dot{Q}-kA\cdot\Delta\vartheta_\mathrm{log} \label{eq:8}
 \end{equation}

In general, the design parameters ($kA$ in case of the heat exchanger) can be
adjusted using lookup table functions to match the model behavior to measured
data. This is especially useful if components with well known characteristics
should be implemented in a different plant or at different operating conditions.
To further improve the representation of an actual plant, the modular structure
of TESPy easily allows the addition of new equations and characteristic
functions to existing components or even new components.

# Previous Implementations

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

# Acknowledgments

This work is supported by University of Applied Sciences and the Center for
Sustainable Energy Systems in Flensburg. It is part of the open energy modeling
framework (oemof), most known for its software package oemof.solph [@oemof].
Many thanks to all
[contributors](https://github.com/oemof/tespy/graphs/contributors).

Key parts of TESPy require the following scientific software packages: CoolProp
[@CoolProp], NumPy [@NumPy], pandas [@pandas]. Other packages implemented are
tabulate and SciPy [@SciPy].

# References
