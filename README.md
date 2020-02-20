# BiomassMW
(The functions in this repository are already incorporated into the **[COBRA toolbox](https://github.com/opencobra/cobratoolbox)** since 2018. It is recommended to install COBRA toolbox to use them.)

MATLAB functions for determining the chemical formulae and molecular weights of macromolecules in genome-scale metabolic models.  
Please see the following paper for more details:  
Siu H. J. Chan, Jingyi Cai, Lin Wang, Margaret N. Simons-Senftle, Costas D. Maranas (2017) Standardizing biomass reactions and ensuring complete mass balance in genome-scale metabolic models, *Bioinformatics*, 33(22), 3603–3609.
**[Link](https://doi.org/10.1093/bioinformatics/btx453)**

Main functions:  
1. `computeMetFormulae` two functionalities:  
- Compute the chemical formulas of the unknown metabolites given the formulas for a set of known metabolites using a set of reactions.
- Compute the minimum and maximum possible MW of a target metabolite by constraining a minimum level of inconsistency in elemental balance.

Other functions:  
1. `getElementalComposition`  
For converting chemical formulas into a matrix and checking the elemental balance of reactions  
2. `elementalMatrixToFormulae`  
For converting a matrix into chemical formulas  

Test script:
`testComputeMetFormulae.m`

Example scripts:
`example1.m`
`example2.m`
