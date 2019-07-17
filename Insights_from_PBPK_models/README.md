# PBPK_host-microbiome_model

Pharmacokinetic Multi-Compartment Modeling

The multi-compartment pharmacokinetic model of drug metabolism in the mouse to separate host and microbial contributions 
to systemic metabolite exposure. The model contains 8 main compartments (small intestine sections jejunum, duodenum and ileum, large intestine sections cecum, colon, and distal colon, gallbladder, and central compartment). Note that the gallbladder compartment is added to the original model for simulation of biliary excretion and enterohepatic circulation.

In this folder, workflows are provided to expand the model simulations to assess the effect of microbiota on the systemic drug and metabolite levels under conditions of altered host physiology, microbiota drug-metabolizing activity or physico-chemical properties of drugs. 

The code is developed with MatLab2017 SimBiology pipeline and is distributed under the terms of the GNU 
General Public License (please read copyright_and_license and LICENSE files for details. 

There are following workflows: 

pub_workflow_Fig2_drug_serum_kinetics.m
This workflow simulates the microbiome effect on the drug and metabolite levels in large intestine and serum at different large intestine permeability parameters.

pub_workflow_Fig3_drug_AtoPratio.m
This workflow simulates the microbiome effect on the systemic drug metabolite exposure for highly absorbed drugs.
 
pub_workflow_Fig4_enterohepatic_parameters.m
This worflow simulates the interplay between enterohepatic circulation coefficients, microbiome de-glucuronidation coefficient and microbiome contribution to drug metabolism. 

These model simulations is part of the work by Maria Zimmermann-Kogadeeva, Michael Zimmermann, and Andrew L. Goodman. 
