# PMI Calculator

This project calculates Process Mass Intensity (PMI) for a chemical synthesis route using the sequence of synthetic steps and step information.  The interactive `Shiny` app allows for interative process specification followed by PMI calculation.

*Draft version, December 2016*

*Jacob Albrecht, Bristol-Myers Squibb*

## To use: 

On the first tab, enter each molecular transformation, specifying the stoichiometry and direct product of each reagent or intermediate in the sequence. The sequence and branch points will automatically be infered from these relationships. If a reagent is not limiting, enter the stoichiometry of the reagent charge for that step.

On the second tab, enter the reagent details for each of the synthetic intermediates- molecular weight and product yield.  Also enter the mass of all other charges to the process (solvents, process aids, and other consumables) per kilo of step product.

On the thrid tab, there are buttons for starting the Monte Carlo calculation, as well as to download the simulaiton results.  Once the calculation is run, the plotted results will automatically be displayed.
