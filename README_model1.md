# PBPK_host-microbiome_model v1

Pharmacokinetic Multi-Compartment Modeling

Model overview – The multi-compartment pharmacokinetic model of drug metabolism in the mouse contains 7 main compartments (small intestine I-III, cecum, colon, feces, and serum).

Two additional compartments (small_intestine_0 and small_intestine) are used as reservoirs for the initial drug dose. The serum compartment incorporats processes occurring in the liver, kidneys and all other body parts apart from the gastrointestinal (GI) tract.

Exposure to the drug is modelled as an input to the small_intestine_serum compartment of the initial amount of drug equal to D*F, where D is the provided dose, and F is the bioavailability coefficient, and input to the small_intestine_0 compartment of the initial amount of D*(1-F).

Drug propagation through the body is driven by the flow of GI material in different GI tract sections and tissue:serum diffusion coefficients.

Model parameters and equations are provided in files example_basic_model_BRV_mut_wt.csv, example_basic_model_BRV_GF_CV.csv, example_basic_model_SRV_GF_CV.csv, example_basic_model_CLZ_GF_CV.csv. The models differ by drug bioavailability and volume of distribution, initiel parameter values and volumes of GI tract tissues between GF (germfree), CV (conventional) and wt and mut monocolonized mice (B. thetaiotaomicron).

Example workflow is provided in pbpk_host_microbiome_workflow_v1. 

Input files are the model and experimental data files (for example, example_basic_model_BRV_mut_wt.csv, example_data_BRV_MUT.csv and example_data_BRV_WT.csv). 

Output files are created (with preffix 'out_' to the input file name) with model fit and predicted values.

The model was created using the MatLab 2017b SimBiology Toolbox (MathWorks).
