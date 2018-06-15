# PMI Calculator

This project calculates Process Mass Intensity (PMI) for a chemical synthesis route using the sequence of synthetic steps and step information.  The interactive [Shiny](http://shiny.rstudio.com/) app allows for interative process specification and PMI calculation.  The app is hosted here: http://acsgcirt-pmi_calculator.shinyapps.io

*June 2018, ACS Green Chemistry Institute Pharmaceutical Roundtable*

## What is Process Mass Intensity?

Process Mass Intensity (PMI) is a measure of a manufactruing process efficiency.  To make a desired chemical compound there are often several potential synthesis routes, but each route will have a different environmental impact.  One way to measure the possible environmental impact of a process is to consider the quantity of material that must be utilized to make the desired product.  PMI is this simple efficiency measure of the ratio of input material to final product:

\\[PMI = \frac{\sum{MassOfMaterials}}{MassOfIsolatedProduct}\\]

Here \\(MassOfMaterials\\) includes the mass of process solvents, chemical reagents, and any other single-use consumables utilized in the process execution.  While PMI cannot distinguish the impact of the individual reagents and consumables used, it provides a simple and accessible means of discriminating among a set of potential processes.

When a chemist considers different retrosynthetic routes they must depend on their experience to know which is likely to be the most efficient.  To avoid personal biases, an analysis of 1800 reactions organized by reaction type and subtype was performed to characterize the range of PMIs reported for specific chemical reaction types.  By defining a sequence of reactions and their corresponding reaction type it is possible to estimate a plausible PMI for any chemical synthesis route using historical information.

# How to use this app: 

For a manufacturing process, the app requests data on two tabs, and then presents plots of the results on the last tab.  Because information is passsed between tabs only when they are clicked, it is best to click through the tabs in order, starting with the "Define Process" tab.

On the "Define Process" tab, enter each molecular transformation, specifying the stoichiometry and direct product of each reagent or intermediate in the process. The sequence and branch points will automatically be infered from these relationships. If a reagent is not limiting, enter the stoichiometry of the reagent charge for that step.  For labeling compounds, the letters A-Z are offered as selections.  However, custom names may be entered for each of the species.  Each transformation is listed as a seperate row, one for each of the tracked inputs.  The process data may be downloaded for reference using the "Download Process Info as RDS file" button.  Similarly, process information may be restored using the upload button.  No user-entered information is logged.

On the second "PMI Values" tab, enter the molecular weight and reaction details for each of the synthetic steps.  A preset dropdown is available to use literature ranges for specific classes of reactions.  For entering custom ranges product yield and step PMI values may be entered. Step PMI is defined as the mass of all inputs to the process (solvents, process aids, and other consumables) per mass of step product.  A graph of the process sequence will be generated to confirm the information entered on the "Define Process" tab.

On the third "Results" tab, there are buttons for starting the Monte Carlo calculation of PMI, as well as the option to download the simulation results.  The checkbox for advanced options allows access for adjusting two simulation options:
 
 - **Number of Monte Carlo simulations** 5000 iterations are used to balance between calculation time and result quality.

 - **Correlation between PMI and Yield** PMI and Yield have a moderate negative correlation from our analysis.  This is used to draw better MC samples.

Once the calculation has been run, the plotted results will automatically be displayed.  A button to download the raw Monte Carlo samples is available if a more customized analysis is needed.  There are three tabs to show different views of the data:

 - **Overall PMI** A histogram of the cumulative PMI for the whole process, 95% intervals are displayed

 - **Step Metrics** A view of the different contributions to PMI.  Distributions of the mass to produce final product are displayed for each of the intermeidates along with the range of total material needed.  

 - **Step Yield vs Step PMI** Plots are generated to show the correlation between PMI and Yield for each of the reaction products.  

*app created by Jacob Albrecht, Bristol-Myers Squibb*

*app concept and design input from Alina Borovika, Jun Li, Martin Eastgate, Bristol-Myers Squibb*