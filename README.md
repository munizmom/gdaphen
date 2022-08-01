# gdaphen

[](https://img.shields.io/badge/Gdaphen-v0.0-pink)
Gdaphen is a R pipeline to identify the most important predictors qualitative and quantitative variables for discrimination of your variable of interest (like genotype or sex or response to treatment), by using phenotypic, ddPCR or clinical data from different diseases. Gdaphen is able to perform the analyses in individual variables and grouped variables that facilitate the assesment of behavioral or diagnostic tests as many variables are scored in the same tests.

This pipeline was developped in [Yann's Herault team](https://www.igbmc.fr/igbmc/missions/annuaire/yann-herault) as part of one of my postDoc projects and thanks to Yann's supervision. Gdaphen came to exist to perform an in deep statistically assessment of scored variables contributing to genotype effects, focusing specifically in elucidating which variables where contributing in explaining gene dosage effects in mouse and rat models carrying genomic duplications or Knock-out models.

Gdaphen makes the analysis of phenotypic data much easier for medical researchers or behaviorists providing an integrated framework to perform: i) pre-processing steps as data imputation or anonymization; ii)  a full statistical assessment based in GLM and RF classifiers implementing a powerful variable selecting method based on MFA to identify which behavioral predictor variables are most important to discriminate your variable of interest; iii) state of the art visualizations ready for publication to support the conclusions of the analyses.

Other contributors that **must** be acknowledged for providing great discussions in stats, in their analytical needs, expectatives and requirements for the plots are **Yann Herault, Claire Gaveriaux-Ruff, Yaping Xue, Celeste Chidiac, Monsier Arnaud Duchon, Gopal Krishna and Valerie Herault**. And the last is not the less important! as Valerie is currently developping further functionalities for this project. 

Thanks also to **Corentin Guioullier** that did his Master in statistics stage with us and tried out the gdaphen scripts in Jupyter notebook.

![gdaphen](https://github.com/munizmom/gdaphen/blob/master/images/graphicalAbstract.jpg)

## Instalation:
```
# Install devtools to make available the function install_github
install.packages("devtools");

#Load the devtools
library(devtools); 

#Install gdaphen
install_github("munizmom/gdaphen");

#Load gdaphen
library(gdaphen);
```

## Packages Dependencies:
* caret
* cowplot
* data.table
* dplyr
* extrafont
* FactoMineR
* gdata
* ggforce
* ggplot2
* ggpubr
* grDevices
* grid
* gridBase
* gridGraphics
* gtable
* gtools
* Hmisc
* lattice
* mlbench
* openxlsx
* PCAmixdata
* randomcoloR
* randomForest
* RColorBrewer
* readxl
* rlist
* Rmisc
* scatterplot3d
* stringr
* tidyr
* xlsx

## gdaphen modules:

Gdaphen provides an environment to perform:

1.  Pre-processing functions to aid in shaping the input data . Stored in the **module preProcessing**. 
2.  Statistical analyses. Stored in the **module analysis**.
3.  Visualizations or plotting functions. Stored in the **module visualization**.

- Note: Please take also a look at the examples provided. These examples will further guide you along the pipeline specially in pre-processing steps, like assigning specific colours to the variables or observations, and will give you examples of how to impute missing values.


## Closer look into gdaphen functions

![Gdaphen functions](https://github.com/munizmom/gdaphen/blob/master/images/gdaphen_functions_1.jpg)

![Gdaphen functions](https://github.com/munizmom/gdaphen/blob/master/images/gdaphen_functions_2.jpg)

![Gdaphen functions](https://github.com/munizmom/gdaphen/blob/master/images/gdaphen_functions_3.jpg)

## Publications using gdaphen:
For now only our team Herault members @IGBMC Strasbourg used it but we hope others... as the reader, will like to implement it.
The gdaphen methods manuscript is undergoing but a detailed explanation of the method can be found [here](https://www.frontiersin.org/articles/10.3389/fphar.2021.780132/full#supplementary-material).


1. [**Chidiac C, Xue Y, Muniz Moreno MDM, et al.** The Human SCN10AG1662S Point Mutation Established in Mice Impacts on Mechanical, Heat, and Cool Sensitivity. *Front Pharmacol*. 2021 Dec 1;12:780132. ](https://www.frontiersin.org/articles/10.3389/fphar.2021.780132/full). [Take a look at figure 5](https://www.frontiersin.org/files/Articles/780132/fphar-12-780132-HTML/image_m/fphar-12-780132-g005.jpg)
2. [**Xue Y, Kremer M, Muniz Moreno MDM, Chidiac C, Lorentz R, Birling MC, Barrot M, Herault Y, Gaveriaux-Ruff C.** The Human SCN9AR185H Point Mutation Induces Pain Hypersensitivity and Spontaneous Pain in Mice. *Front Mol Neurosci.* 2022 Jun 13;15:913990. PMID: 35769334.](https://www.frontiersin.org/articles/10.3389/fnmol.2022.913990/full)
3. [**Duchon A, Muniz Moreno MdM, Chevalier C, Nalesso V, Andre P, Fructuoso-Castellar M, Mondino M, Po C, Noblet V, Birling MC, Potier MC and Herault Y.** www.biorxiv.org](https://www.biorxiv.org/content/10.1101/2022.06.06.494940v1.full.pdf)


## Contributions
Pull requests are always welcomed!. For major changes, problems, extra functionalities or if you need some help to run gdaphen in your own data please open an issue first or contact  Yann Herault (herault@igbmc.fr) and me (munizmorenomariadelmar@gmail.com) to discuss.

## License
[GPL-3](https://www.gnu.org/licenses/gpl-3.0.html)


