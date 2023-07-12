ECOTRAN is an end-to-end ecosystem modelling platform developed by Steele and Ruzicka (2011) and further developed by Jim Ruzicka and others (see reference list at end). 

README_ECOTRAN_NCC2_05182023.docx is a more comprehensive readme for the ECOTRAN modelling platform and the functions existing within this repository.

"ECOTRANdynamic_NCC2_TEMPLATE_PARALLEL_05152023.m" is the only file the user needs to open if they want to run scenarios with the food web model "NCC2_09032022.csv". 

Open "ECOTRANdynamic_NCC2_TEMPLATE_PARALLEL_05152023.m" and the first lines of code can be edited to run different scenarios. These lines are annotated to explain each section to the user. When all lines of code are edited, run this script, and all other functions will be called internally. Direct any questions to Dylan Gomes (dylan.ge.gomes@gmail.com) or Jim Ruzicka (James.Ruzicka@noaa.gov).


For more advanced uses, "NCC2_09032022.xlsx" can be modified if the user would like to change functional group parameters. Within this, there is a macro button on the "Main" tab that exports the file to the "NCC2_09032022.csv" file. Excel export locations can vary by user settings (defaults often to "Documents/").

"NCC_SubRegional/" shows files for breaking out the subregions, but this isn't necessary to run scenarios. "ECOTRANdynamic_NCC2_TEMPLATE_PARALLEL_05152023.m" can allow the user to select subregions within the first few lines of code with the `Region` object.




References:

Steele, John H., and James J. Ruzicka. "Constructing end-to-end models using ECOPATH data." Journal of Marine Systems 87.3-4 (2011): 227-238. https://doi.org/10.1016/j.jmarsys.2011.04.005

Ruzicka, J. J., Brink, K. H., Gifford, D. J., & Bahr, F. (2016). A physically coupled end-to-end model platform for coastal ecosystems: Simulating the effects of climate change and changing upwelling characteristics on the Northern California Current ecosystem. Ecological Modelling, 331, 86–99. https://doi.org/10.1016/j.ecolmodel.2016.01.018

Ruzicka, J. J., Brodeur, R. D., Emmett, R. L., Steele, J. H., Zamon, J. E., Morgan, C. A., Thomas, A. C., & Wainwright, T. C. (2012). Interannual variability in the Northern California Current food web structure: changes in energy flow pathways and the role of forage fish, euphausiids, and jellyfish. Progress in Oceanography, 102, 19-41. https://doi.org/10.1016/j.pocean.2012.02.002

Ruzicka, J. J., Daly, E. A., & Brodeur, R. D. (2016). Evidence that summer jellyfish blooms impact Pacific Northwest salmon production. Ecosphere, 7(4). https://doi.org/10.1002/ecs2.1324 

Ruzicka, J. J., Kasperski, S., Zador, S., & Himes-Cornell, A. (2019). Comparing the roles of Pacific halibut and arrowtooth flounder within the Gulf of Alaska ecosystem and fishing economy. Fisheries Oceanography, 28, 576-596. https://doi.org/10.1111/fog.12431

Ruzicka, J. J., Steele, J. H., Ballerini, T., Gaichas, S. K., & Ainley, D. G. (2013). Dividing up the pie: whales, fish, and humans as competitors. Progress in Oceanography, 116, 207-219. https://doi.org/10.1016/j.pocean.2013.07.009

Ruzicka, J. J., Steele, J. H., Gaichas, S. K., Ballerini, T., Gifford, D. J., Brodeur, R. D., & Hofmann, E. E. (2013). Analysis of energy flow in US GLOBEC ecosystems using end-to-end models. Oceanography, 26(4), 82-97. http://dx.doi.org/10.5670/oceanog.2013.77. 

Ruzicka, J. J., Steele, J. S., Brink, K. H., Gifford, D. J., & Bahr, F. (2018). Understanding large-scale energy flows through end-to-end shelf ecosystems - the importance of physical context. Journal of Marine Systems, 187, 235-249. https://doi.org/10.1016/j.jmarsys.2018.08.003

Gomes, D., Ruzicka, J. J., Crozier, L. G., Huff, D. D., Phillips, E. M., Hernvann, P. Y., ... & Auth, T. D. (2022). An updated end-to-end ecosystem model of the Northern California Current reflecting ecosystem changes due to recent marine heat waves. bioRxiv, 2022-12. https://doi.org/10.1101/2022.12.28.522165

