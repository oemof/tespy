Benchmarks
==========
To ensure credibility as well as reproducibility for the software, several
measures are taken. The most important information is listed below:

Model Validation
----------------
TESPy has been used to model several research and engineering applications. In
the paper on integration of generic exergy analysis in TESPy
:cite:`Witte2022` three models have been built from literature sources: A
solar thermal power plant, a supercritical CO2 Brayton cycle as well as a
refrigeration machine using air as working fluid.

For the solar thermal power plant we have created a full model of the plant
using a standard industry software in parallel. **The comparison showed**
**identical results**. For the other two applications we have compared the
results of the TESPy model with the data published in the respective research
paper and found very well matching results. Differences can be explained by
different implementations of the fluid property back-end.

Finally, in the extension of the exergy analysis to chemical exergy
:cite:`Hofmann2022` we have also compared results of the CGAM process
:cite:`Valero1994` modeled in TESPy with a full model using industry software
and with the data provided from literature as well :cite:`Bejan1996`.

The code for the full models is accessible open source on GitHub:

- `So called "Solar Energy Generating System" <https://github.com/fwitte/SEGS_exergy>`__
- `Supercritical CO2 power cycle <https://github.com/fwitte/sCO2_exergy>`__
- `Air refrigeration machine <https://github.com/fwitte/refrigeration_cycle_exergy>`__
- `CGAM process <https://github.com/KarimHShawky/Chemical-Exergy-in-TESPy>`__


Unit testing
------------
On top of the full model validation, the software includes full unit testing.
Here, single features of the software are tested by comparing the result the
software provides with the result we would expect when manually modeling that
feature. For example, we set up a turbine and check, whether the isentropic
efficiency specification in the TESPy component matches the results that is
expected doing the same process manually. This is done for all modules of the
software.

Continuous Integration
----------------------
TESPy has a
`Continuous Integration <https://en.wikipedia.org/wiki/Continuous_integration>`__
pipeline. The unit tests and the full model tests are automatically run
whenever changes to the source code of the software are made. By this we can
ensure, that changes in the code do not break with existing features or
invalidate results of the existing models.

For more information on how to run the tests please see the
:ref:`how to develop <tespy_development_how_label>` section.
