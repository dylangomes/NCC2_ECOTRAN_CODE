"EcoTran_Code/" contains all the files necessary to run simulations. "ECOTRANdynamic_NCC2_TEMPLATE_PARALLEL_05152023.m" is the only file the user needs to open if they want to run scenarios with the food web model "NCC2_09032022.csv". 

Open "ECOTRANdynamic_NCC2_TEMPLATE_PARALLEL_05152023.m" and the first lines of code can be edited to run different scenarios. These lines are annotated to explain each section to the user. When all lines of code are edited, run this script, and all other functions will be called internally.



For more advanced uses, "NCC2_09032022.xlsx" can be modified if the user would like to change functional group parameters. Within this, there is a macro button on the "Main" tab that exports the file to the "NCC2_09032022.csv" file. Excel export locations can vary by user settings (defaults often to "Documents/").

"NCC_SubRegional/" shows files for breaking out the subregions, but this isn't necessary to run scenarios. "EcoTran_Code/ECOTRAN_DynamicScenarios_NCC2_2D_Par.m" can allow the user to select subregions in the first few lines of code.

See "README_ECOTRAN-2_05132021.doc" for a description of ECOTRAN model code