
#' mfa.function - Function to perform Multiple Factor Analysis for mixed data (MFAmix) 
#'
#' @name mfa.function 
#'
#'
#' @param dataInputPCA Datframe where each row represent an animal/pair of animal measurements, and the columns represents  
#'  qualitative/quantitative parameters measurements. Data should be complete no NAs accepted, 
#'  if needed imputation strategies should be used before running this function.
#'
#' @param catVarNb Number of categorical variables contained in dataInputPCA. Ex use 2 if Genotype and Sex  are present.
#'
#' @param nameOutput Sufix Name specific from this project to add to the name of the output of the function
#'
#' @return A list of dataframes object: [res_PCAmixData_]nameOutput

#' @examples
#' # for discrimination of  Genotype:
#' mfa.function(parameters_more_30perc_informative_Genotype,2, "parameters_more_30perc_informative_Genotype");
#' # after removing Genotype, which variables are the most discriminative?
#' mfa.function(parameters_more_30perc_informative_woGenotype_Genotype,1, "parameters_more_30perc_informative_woGenotype");
#' # for discrimination of  Sex:
#' mfa.function(parameters_more_30perc_informative_Sex,2, "parameters_more_30perc_informative_Sex");
#'
#' @export
mfa.function <- function(dataInputPCA,catVarNb, nameOutput){
	# description: Function to perform Multiple Factor Analysis for mixed data (qualitative and quantitative variables) (using MFAmix) 
	#library("PCAmix data");
	########################################################
	# results of step 4: MFA to feed to the plot functions
	########################################################
	#-> output file: res_PCAmixData_phenoAll_4analysis
	##  Results of the Multiple Factor Analysis for mixed data (MFAmix) **

	## The analysis was performed on x individuals, described by xx variables
	## Results are available in the following objects :
	####################################################
	#"name" "description"
	#"$eig" "eigenvalues"
	#"$eig.separate" "eigenvalues of the separate analyses"
	#"$separate.analyses" "separate analyses for each group of variables"
	#"$groups" "results for all the groups"
	#"$partial.axes" "results for the partial axes"
	#"$ind" "results for the individuals"
	#"$ind.partial" "results for the partial individuals"
	#"$quanti" "results for the quantitative variables"
	#"$levels" "results for the levels of the qualitative variables"
	#"$quali" "results for the qualitative variables"
	#"$sqload" "squared loadings"
	#"$listvar.group" "list of variables in each group"
	#"$global.pca" "results for the global PCA"

	names <- unique(as.vector(unlist(gsub("\\:: .*","", colnames(dataInputPCA)))));
	nbPhenos <- data.frame(phenoEvents=unique(names));
	nbPhenos$class.var <- 1:nrow(nbPhenos);
	namesDf <- data.frame(phenoEvents=as.vector(unlist(gsub("\\:: .*","", colnames(dataInputPCA)))));
	namesDf <- left_join(namesDf,nbPhenos,by="phenoEvents");
	class.var <- as.vector(unlist(namesDf$class.var));
	

	maxNb_FactorsOutcomes.tmp <- as_tibble(data.frame(terms=as.vector(unlist(gsub("\\:: .*","", colnames(dataInputPCA))))));
	maxNb_FactorsOutcomes <- max(as.data.frame(maxNb_FactorsOutcomes.tmp %>% dplyr::group_by(terms) %>% dplyr::summarise(n = n()))[,2])
	
	if (maxNb_FactorsOutcomes < 3){
		maxNb_FactorsOutcomes <- 3;
	};

	pdf(paste0( wd,"/",f_pheno_discrimination,f_results,nameOutput,"_PCAmixData.pdf"));
	res<-MFAmix(data=dataInputPCA,groups=class.var,name.groups=names, rename.level=TRUE, ndim=maxNb_FactorsOutcomes,graph=TRUE);
	summary(res);
	dev.off();

	assign(paste0("res_PCAmixData_",nameOutput),res,.GlobalEnv);

}; 


#' More informativeVars in PCA - Identification of Sel30 variables
#'
#' @name informative_params_more_30.function 
#'
#' @param InputdataPCAMix The input is the output from the mfa.function, the object: [res_PCAmixData_]nameOutput 
#'
#' @param catVarName Column name of the qualitative
#'  variable or categorical variable that is your dependant 
#'  variable, the one for which one you want to identify the most discriminative variables  
#'
#' @param phenoColourDf Dataframe containing at least 2 columns, 
#' 1: "labels", containing the possible categories. Ex Genotype het, wt, homo
#' 2: "colPalette" hex colour associated to each varible. 
#' Ex: for label="heterozygous", use colPalette="#8B0000" 
#' 
#' @return A dataframe, like dataInputPCA but containing only the columns
#' or variables thata re contributing more than a 30% to the discrimination
#' Note: Sex and genotype are always preselected.These columns do not follow
#' the order the experimentator provided in phenoSocial_4analysis. The output has a prefix
#' catVarName, and followed by "_sel_Vars_contribution_more_30perc_plus_categorical"
#'
#' @examples
#' # for discrimination of  Genotype
#' informative_params_more_30.function(res_PCAmixData_phenoAll_4analysis_woHighCor_sel1,"Genotype",phenoColourDf);
#' # for discrimination of  Sex
#' informative_params_more_30.function(res_PCAmixData_phenoAll_4analysis_woHighCor_sel1,"Sex",phenoColourDf);
#'
#' @export
informative_params_more_30.function <- function(InputdataPCAMix,catVarName,phenoColourDf){
	# description: For your variable of interest, selection of explanatory features 
	# contributing to the discrimination more than a 30% after running the MFA 
	# analysis.
	#
    #df1: list to df
	for (i in 1:length(InputdataPCAMix$partial.axes)){
		if (i==1){
			dfdata1 <- as.data.frame(InputdataPCAMix$partial.axes[[1]]);
			df1 <- dfdata1;
		} else {
			dfdata.tmp <- as.data.frame(InputdataPCAMix$partial.axes[[i]]);
			df1 <- rbind(df1, dfdata.tmp);
			};
			i=i+1;
	};


	labels_full <- rownames(df1); #plot1: Grouped-variables correlations  
	genotypesNames <- gsub("^.*=","",rownames(InputdataPCAMix$levels$coord)[grep(catVarName, rownames(InputdataPCAMix$levels$coord))])
	genotype_last_dimension <- length(labels_full[grep(catVarName,labels_full)])
	genotypesNames <- genotypesNames[c(1:genotype_last_dimension)]	    
	Df1_wogenotypesLabels <- labels_full[-grep(catVarName,labels_full)]
    Df1genotypesLabels <- paste0(labels_full[grep(catVarName,labels_full)],"_",genotypesNames)
    labels_full <- c(Df1genotypesLabels,Df1_wogenotypesLabels )
    rownames(df1) <- labels_full; #plot1: Grouped-variables correlations  
    df1$labels <- gsub(".*\\.","",rownames(df1)); #plot1: Grouped-variables correlations  

	df1<- left_join(df1,phenoColourDf,by="labels" );
	df1$nb <- 1:nrow(df1);
	df1$labels_nb <- paste0(df1$nb,") ",df1$labels);
    df1$labels_full <- labels_full; #plot1: Grouped-variables correlations  

	#all parameters tested and order
	if (length(genotypesNames)==1){
		dfnames <- as.data.frame(unlist(InputdataPCAMix$listvar.group));
		namesRows <- rownames(dfnames);	
		colnames(dfnames) <- "parameters_more_30perc_informative";
		dfnames$nb <- as.numeric(1:nrow(dfnames));

	} else{
	dfnames <- as.data.frame(unlist(InputdataPCAMix$listvar.group));
	dfnames[(nrow(dfnames)+(length(genotypesNames)-1)),] <-dfnames[1,];
	namesRows <- rownames(dfnames);	
	namesRows<- 	namesRows[c(length(namesRows), 1:(length(namesRows)-1))];
	namesRows[c(1:length(genotypesNames))] <- paste0(catVarName,"_",genotypesNames);
	dfnames<- 	as.data.frame(dfnames[c(nrow(dfnames), 1:(nrow(dfnames)-1)),]);
	rownames(dfnames)<- namesRows;
	colnames(dfnames) <- "parameters_more_30perc_informative";
	dfnames$nb <- as.numeric(1:nrow(dfnames));
	}
	
 	#Second selection: most informative variables, selecton to run the analysis again based on the degree of info 
	# Lets select now the most informative variables to do another round of selections, by setting up the threshold of selection
	# only the variables contributing in the discrimination more than 0.3 or a 30% will be used plus sex and genotype than 
	# in no case we want to lose
    #Creating the input files for tle 6 plots:
    if ("Sex" %in% df1[,'labels']==TRUE){
		sel_Vars_contribution_more_30perc_plus_categorical <- unique(rbind(df1[grep(catVarName, df1[,'labels']),],
			df1[which(df1[,'labels']=="Sex"),],df1[which(df1[,1]>0.3),],df1[which(df1[,2]>0.3),],df1[which(df1[,3]>0.3),]));
	} else {
		sel_Vars_contribution_more_30perc_plus_categorical <- unique(rbind(df1[grep(catVarName,df1[,'labels']),],
			df1[which(df1[,1]>0.3),],df1[which(df1[,2]>0.3),],df1[which(df1[,3]>0.3),]));
	};
	sel_Vars_contribution_more_30perc_plus_categorical <- left_join(sel_Vars_contribution_more_30perc_plus_categorical,dfnames,by="nb")

    #assign("dfnames",dfnames,.GlobalEnv);
    assign(paste0(catVarName,"_sel_Vars_contribution_more_30perc_plus_categorical"),sel_Vars_contribution_more_30perc_plus_categorical,.GlobalEnv);
	
};



#' Parameters_more_30perc_informative - Creating the Sel30 model
#'
#' @name parameters_more_30perc_informative_creation_InputDf.function 
#'
#' @param phenoSocial_4analysis The input is the output from the mfa.function, the object: [res_PCAmixData_]nameOutput 
#'
#' @param sel_Vars_contribution_more_30perc_plus_categorical Output of the function
#' informative_params_more_30.function
#'  variable or categorical variable that is your dependant 
#'  variable, the one for which one you want to identify the most discriminative variables  
#'
#' @param nameOutput Sufix Name specific from this project to add to the name of the output of the function
#' 
#' @return A dataframe, like phenoSocial_4analysis but containing only the columns
#' or variables that are contributing more than a 30% to the discrimination
#' Contrary to the dataframe produced by 
#' informative_params_more_30.function
#' in this case, the columns are following the order specified by the experimentator in phenoSocial_4analysis
#' These columns do not follow the order the experimentator provided in the dataInputPCA. 
#' The output has  aprefix= "parameters_more_30perc_informative_" and as  a sufix nameOutput
#'
#' @examples
#' # for discrimination of  Genotype:
#' parameters_more_30perc_informative_creation_InputDf.function(phenoSocial_4analysis,Genotype_sel_Vars_contribution_more_30perc_plus_categorical,"Genotype");
#' # for discrimination of  Sex:
#' parameters_more_30perc_informative_creation_InputDf.function(phenoSocial_4analysis,Sex_sel_Vars_contribution_more_30perc_plus_categorical,"Sex");
#'
#' @export
parameters_more_30perc_informative_creation_InputDf.function<- function(phenoSocial_4analysis,sel_Vars_contribution_more_30perc_plus_categorical,nameOutput){
	# description: For your variable of interest, selection of explanatory features 
	# contributing to the discrimination more than a 30% after running the MFA 
	# analysis.
	
	parameters_more_30perc_informative <-phenoSocial_4analysis[, (colnames(phenoSocial_4analysis) %in% sel_Vars_contribution_more_30perc_plus_categorical$parameters_more_30perc_informative)]; 
	parameters_more_30perc_informative_woGenotype<- parameters_more_30perc_informative[ , !(names(parameters_more_30perc_informative) %in% "Genotype")];
	assign(paste0("parameters_more_30perc_informative_",nameOutput),parameters_more_30perc_informative,.GlobalEnv);
	assign(paste0("parameters_more_30perc_informative_woGenotype_",nameOutput),parameters_more_30perc_informative_woGenotype,.GlobalEnv);

};


#' MFA extraction - Creation of plots input dataframes 
#'
#' @name any_categorical_variable_discrimination_byPhenoParams.function 
#'
#' @param InputdataPCAMix The input is the output from the mfa.function, the object: 
#'[res_PCAmixData_]nameOutput 
#'
#' @param nameExperiment The name of the experiment together with the number of variables we included. 
#' if is all model, allthe variables. Instead if is the selected model is less variables. 
#' ex: 
#' 1. nameExperiment <- paste0(name,"_", length(colnames(phenoSocial_4analysis)), " variables"); #all variables
#' 2. nameExperiment <- paste0(name,"_", length(colnames(phenoAll_4analysis_woHighCor_sel1)), " variables"); · without highly correlated
#' 3. nameExperiment <- paste0(name,"_", " more informative variables ",length(colnames(parameters_more_30perc_informative_Genotype)) ) #selected variables
#'
#' @param catVarName Column name of the qualitative
#'  variable or categorical variable that is your dependant 
#'  variable, the one for which one you want to identify the most discriminative variables  
#' 
#' @param phenoColourDf Dataframe containing at least 2 columns, 
#' 1: "labels", containing the possible categories. Ex Genotype het, wt, homo
#' 2: "colPalette" hex colour associated to each varible. 
#' Ex: for label="heterozygous", use colPalette="#8B0000" 
#' 
#' @param indSexCond Samples data color asignation for each Individual point depending
#' on genotypes, sex, treatments is created by the user.
#' an example of how to create it is in the auxiliar function
#' ' indSexCond.colour_definition.function.
#'  
#' @param InputDf Is always the phenoSocial_4analysis file created in the preproccesing steps 
#' condition_definition.function. 
#' 
#' @param qualiVars Number of qualitative variables used on the model.
#'
#' @param nameModel Shrt naming for the model you are using.
#' Ex: "fullMdl" for the model including all variables;"Sel1" for the model not including the
#' highly correlated variables, "Sel30" for the model including only the variables contributing
#' more thana 30%
#' 
#' @return 6 dataframe and several variables that contain the information for plotting
#' namesVariablesDf, df1, df2, df3,df4, df5, df6, MFA1, MFA2, MFA3,nbVarF, nbVarTotalF'
#'
#' @examples
#'any_categorical_variable_discrimination_byPhenoParams.function(InputdataPCAMix,nameExperiment,catVarName,phenoColourDf,indSexCond,phenoSocial_4analysis,qualiVars,nameModel)
#' 
#' @export
any_categorical_variable_discrimination_byPhenoParams.function <- function(InputdataPCAMix,nameExperiment,catVarName,phenoColourDf,indSexCond,InputDf,qualiVars,nameModel){
 	# description: Formatting the mfa.function output (InputdataPCAMix).
	# 1. Storing the info of the mfa.function in input dataframes, one for each plot we will compute. 
	# 2. Then formating those dataframes (df1, for plot1, df2 for plot2, etc 
	# to include plot labels and the specific colours to use for the plotting.
	#
   #df1: list to df
	for (i in 1:length(InputdataPCAMix$partial.axes)){
		if (i==1){
			dfdata1 <- as.data.frame(InputdataPCAMix$partial.axes[[1]]);
			df1 <- dfdata1;
		} else {
			dfdata.tmp <- as.data.frame(InputdataPCAMix$partial.axes[[i]]);
			df1 <- rbind(df1, dfdata.tmp);
			};
			i=i+1;
	};

    labels_full <- rownames(df1); #plot1: Grouped-variables correlations  
	print(labels_full);
	labelsQuanti<- gsub(':(?!.*:) ', ': (', gsub(" $",")",gsub("\\."," ",gsub("\\.\\.",": ",gsub("\\.\\.\\.",":: ",rownames(InputdataPCAMix$quanti$coord))))),perl = TRUE);

	for (i in 1:length(qualiVars)){
		if (i==1){
			qualiVarsNb <- length(grep(qualiVars[i], labels_full));
		} else {
			qualiVarsNb.tmp <- length(grep(qualiVars[i], labels_full));
			qualiVarsNb <- qualiVarsNb +qualiVarsNb.tmp;
		};
		i= i+1;
	};

	beginingQuantiVars <- length(labels_full)-(length(labels_full)- qualiVarsNb) +1
	dimInfo <- gsub("\\..*$","",labels_full[beginingQuantiVars: length(labels_full)]); 
	labels_full[beginingQuantiVars: length(labels_full)] <- paste0(dimInfo,".",labelsQuanti)

	genotypesNames <- gsub("^.*=","",rownames(InputdataPCAMix$levels$coord)[grep(catVarName, rownames(InputdataPCAMix$levels$coord))])
	genotype_last_dimension <- length(labels_full[grep(catVarName,labels_full)])
	genotypesNames <- genotypesNames[c(1:genotype_last_dimension)]	    

    Df1genotypesLabels <- paste0(labels_full[grep(catVarName,labels_full)],"_",genotypesNames)
	labels_full[1: length(Df1genotypesLabels)] <- Df1genotypesLabels
    rownames(df1) <- labels_full; #plot1: Grouped-variables correlations  

	df1$nb <- 1:nrow(df1);
    df1$labels_full <- labels_full; #plot1: Grouped-variables correlations  
	df1$labels_nb <- paste0(1:nrow(df1), ") ",gsub("_.\\.","_",df1$labels_full));
	df1$labels_nb <- gsub("dim .*\\.","",df1$labels_nb);
    df1$labels <- gsub("::.*","",gsub("Genotype_","",gsub("^.*\\) ","",df1$labels_nb))); #plot1: Grouped-variables correlations  
	df1<- left_join(df1,phenoColourDf,by="labels" );

	namesVariablesDf <- data.frame(Variable=colnames(InputDf));
	namesVariablesDf$labels <- gsub("::.*","",namesVariablesDf$Variable);
	namesVariablesDf<- left_join(phenoColourDf,namesVariablesDf,by="labels") 
	namesVariablesDf[is.na(namesVariablesDf$Variable),'Variable'] <- namesVariablesDf[is.na(namesVariablesDf$Variable),'labels']

	df2 <- as.data.frame(InputdataPCAMix$groups$contrib); #plot2:Grouped variables contribution/correlation to the 3 main PCAs 
	df2$labels <- gsub("\\."," ",gsub("Social.","Social:",gsub("\\.\\.\\..*","",rownames(df2))));
	df2<- left_join(df2,phenoColourDf,by="labels" );

	df3 <- as.data.frame(InputdataPCAMix$ind$coord); # plot3: Individual component map PCA or Individual observations distribution 
	#                                                # to the classifier discimination of qualitative variables in the PCs 
	df3$labels <- rownames(df3);
	df3<- left_join(df3,indSexCond,by="labels" );
	colnames(df3) <- gsub("::.*","",colnames(df3));
	df4 <- as.data.frame(InputdataPCAMix$sqload); ##plot4: Squared loads or Squared loadings of the quantitative & qualitative Vars and the principal components 
	df4$labels <- gsub("\\."," ",gsub("Social.","Social:",gsub("\\.\\.",": ",gsub("\\.\\.\\..*","",rownames(df4)))));
	df4<- left_join(df4,phenoColourDf,by="labels" );

	df5 <- as.data.frame(InputdataPCAMix$quanti$coord);#plot5: Correlation circle  of all variables
	df5$labels <- gsub("\\."," ",gsub("Social.","Social:",gsub("\\.\\.",": ",gsub("\\.\\.\\..*","",rownames(df5)))));
	df5$labels_full <- gsub(" $","",gsub("\\."," ",gsub("::::::","",gsub("\\::\\.$","",gsub("\\.\\.",": ",gsub("\\.\\.\\.","::",rownames(df5)))))));
	df5$labelsParameter <- gsub(" $","",gsub("\\."," ",gsub(".*\\.\\.","",gsub("\\.\\.nb","",
		gsub("\\.\\.g","",gsub("\\.\\.s","",gsub("\\.\\.\\.",": ",gsub("\\.$","",gsub("\\.........$","",rownames(df5))))))))));
	df5$nb <- 1:nrow(df5)
	df5$labels_nb <- paste0(1:nrow(df5), ") ",gsub("_.\\.","_",df5$labels_full)); 
	df5<- left_join(df5,phenoColourDf,by="labels" );

	df6 <- as.data.frame(InputdataPCAMix$levels$coord); #plot6: Levels component map or Qualitative Vars discrimination: Variables distribution in the main dimentional spaces 
	df6$labels <- gsub(".*\\=","",rownames(df6)); 
	df6$Variable <- gsub("\\."," ",gsub("Social.","Social:",gsub("\\.\\.\\..*","",rownames(df6))));
	df6<- left_join(df6,phenoColourDf,by="labels");
	#percentages of variance explained by each of the 3 main components
	MFA1 <-round(InputdataPCAMix$eig[1,'Proportion'] , digits=2);
	MFA2 <-round(InputdataPCAMix$eig[2,'Proportion'] , digits=2);
	MFA3 <-round(InputdataPCAMix$eig[3,'Proportion'] , digits=2); 

	#nbVariables used to build the plots by the PCAMix function nbVarF
	for (i in 1:length(InputdataPCAMix$listvar.group)){
		if (i==1){
			nbVar1 <- as.data.frame(InputdataPCAMix$listvar.group[[1]]);
			nbVar <- nbVar1;
			colnames(nbVar)<- "Variables"
		} else {
			nbVar1.tmp <- as.data.frame(InputdataPCAMix$listvar.group[[i]]);
			colnames(nbVar1.tmp)<- "Variables"
			nbVar <- rbind(nbVar, nbVar1.tmp);
		};
		i=i+1;
	};

	nbVarF<- nrow(nbVar)
	#nb of total variables added to the analyses
	nbVarTotalF <- ncol(InputDf)

    save(namesVariablesDf, df1, df2, df3,df4, df5, df6, MFA1, MFA2, MFA3,nbVarF, nbVarTotalF, file=paste0(wd, f_pheno_discrimination,f_results, f_Rdata,name,"_",nameModel,"_",catVarName,"_Robject_for_plots.RData"))
	save(df3, file=paste0(wd, f_pheno_discrimination,f_results, f_Rdata,name,"_",nameModel,"_",catVarName,"_df3_Robject_for_plots.RData"))
};



#' Linear equation 
#'
#' @name linear_equation_phenoVars_discrimination.function 
#'
#' @param InputdataPCAMix The input is the output from the mfa.function, the object: 
#'[res_PCAmixData_]nameOutput 
#'
#' @param nameExperiment The name of the experiment together with the number of variables we included. 
#' if is all model, allthe variables. Instead if is the selected model is less variables. 
#' ex: 
#' 1. nameExperiment <- paste0(name,"_", length(colnames(phenoSocial_4analysis)), " variables"); #all variables
#' 2. nameExperiment <- paste0(name,"_", length(colnames(phenoAll_4analysis_woHighCor_sel1)), " variables"); · without highly correlated
#' 3. nameExperiment <- paste0(name,"_", " more informative variables ",length(colnames(parameters_more_30perc_informative_Genotype)) ) #selected variables
#'
#' @return A dataframe saved in excel with the linear equation computed.
#'
#' @examples
#' linear_equation_phenoVars_discrimination.function(InputdataPCAMix,nameExperiment);
#' 
#' @export
linear_equation_phenoVars_discrimination.function <- function(InputdataPCAMix,nameExperiment){
	# description: calculate the linear equation vased on the independant variables 
	# considered in the model
	
	resLinear <- data.frame(coefficients=InputdataPCAMix$coef, Variable=rownames(InputdataPCAMix$coef$dim2))
    write.xlsx(resLinear, file = paste0(wd,f_pheno_discrimination, f_results,f_linear_equation,nameExperiment, "_linearEquation_slopes_behavioralVariables.xlsx"));  
};


#' Cosine Distance - Function to calculate cosine simmilarity distances
#'
#' @name stats_cosineDistance_plot1.function 
#'
#' @param nameExperiment The name of the experiment together with the number of variables we included. 
#' if is all model, allthe variables. Instead if is the selected model is less variables. 
#' ex: 
#' 1. nameExperiment <- paste0(name,"_", length(colnames(phenoSocial_4analysis)), " variables"); #all variables
#' 2. nameExperiment <- paste0(name,"_", length(colnames(phenoAll_4analysis_woHighCor_sel1)), " variables"); · without highly correlated
#' 3. nameExperiment <- paste0(name,"_", " more informative variables ",length(colnames(parameters_more_30perc_informative_Genotype)) ) #selected variables
#'
#' @param catVarName Column name of the qualitative
#'  variable or categorical variable that is your dependant 
#'  variable, the one for which one you want to identify the most discriminative variables  
#' 
#' @param similarityThreshold to identify the simmilarity in cosine distance, this parameter 
#' is the cut of simmilarity. Above this threshold we consider the variables highly similar.
#' Ex: 0.75, or 0.85
#' 
#' @return  calculation of cosine simmilarity distances
#'
#' @examples
#' similarityThreshold <- 0.86
#' stats_cosineDistance_plot1.function( catVarName,similarityThreshold, nameExperiment);
#'
#' @export
stats_cosineDistance_plot1.function <- function(catVarName,similarityThreshold, nameExperiment){
	#
	# description: calculation of cosine distances to produce the tables to interpret plot1 Formatting the mfa.function output (InputdataPCAMix).
 	####################### Creating the folders needed ################
	f_similarityThreshold <- paste0(similarityThreshold,"/")
	f_cosineDistance <- "cosineDistance/"
	f_similarityThreshold <- paste0(similarityThreshold,"/")
	dir.create(file.path(wd, f_pheno_discrimination, f_results, f_cosineDistance), showWarnings = F);
	dir.create(file.path(wd, f_pheno_discrimination, f_results, f_cosineDistance,f_similarityThreshold), showWarnings = F);

 	####################### Defining the functions needed ################
	extract_cosSimilarity_vars_using_threshold.function <- function(cosDist, catVarName, similarityThreshold,inputDfNam,name_cosDist,nameExperiment){

		resList_plot<- list();
		j<-1
		for (i in 1:ncol(cosDist)) {
			res <- cosDist[,i][which(cosDist[,i]>= similarityThreshold)];
			if ( length(res)>0){
				res <- as.data.frame(res);
				colnames(res) <-colnames(cosDist)[i];
				resList_plot[[j]] <- res
				names(resList_plot)[j] <- colnames(cosDist)[i];
				j<- j+1;
			};
			i <- i+1;
		};
		if (length(resList_plot)>0){
			write.xlsx(resList_plot, file = paste0(wd, f_pheno_discrimination, f_results,f_cosineDistance,f_similarityThreshold,nameExperiment,"_",catVarName,"_cosSim_",
				name_cosDist,".xlsx"),colNames = TRUE,rowNames = TRUE,);
			assign(paste0(catVarName,"_selectedVars_cosDistance_",name_cosDist,gsub("\\.","-",similarityThreshold)),resList_plot,.GlobalEnv);
			
		};
	};
		####################################################################

 	####################################################################
	# Loading the neccesary libraries and input files
	####################################################################
	#library("SnowballC");
	#library("lsa");

    load( file=paste0(wd, f_pheno_discrimination,f_results, f_Rdata,name,"_",catVarName,"_Robject_for_plots.RData")) 
    statsInputDf<- df1;
    #RData obejct storing: namesVariablesDf, df1, df2, df3,df4, df5, df6, MFA1, MFA2, MFA3,nbVarF, nbVarTotalF,
 	####################################################################
   
 	####################################################################
    # Pre-processing of the input file to reorder cat variables first
	# Stats calculation: cosine distance and evctor lengths correlation
   

	DfStats <- statsInputDf[, colnames(statsInputDf)[grep("dim",colnames(statsInputDf))]];
	rownames(DfStats) <-statsInputDf$labels_full;
	varsName <- statsInputDf$labels_full;
	varsNames_considered_rows<- varsName[-grep(catVarName,varsName)]; # all the other quali and quantitative avriables except the oens coming from the  catVarName we want to know smt about
	varsNames_catVarName<- varsName[grep(catVarName,varsName)]; #genotype effect, line or sex effect

	cosDistanceTotal <- cosine(as.matrix(t(DfStats)));
	colnames(cosDistanceTotal) <- statsInputDf$labels_full;
	rownames(cosDistanceTotal) <- statsInputDf$labels_full;
	
	cosDistance_dim1_vs_dim2 <- cosine(as.matrix(t(DfStats[,c("dim 1", "dim 2")])));
	colnames(cosDistance_dim1_vs_dim2) <- statsInputDf$labels_full;
	rownames(cosDistance_dim1_vs_dim2) <- statsInputDf$labels_full;

	cosDistance_dim1_vs_dim3 <- cosine(as.matrix(t(DfStats[,c("dim 1", "dim 3")])));
	colnames(cosDistance_dim1_vs_dim3) <- statsInputDf$labels_full;
	rownames(cosDistance_dim1_vs_dim3) <- statsInputDf$labels_full;

	cosDistance_dim2_vs_dim3 <- cosine(as.matrix(t(DfStats[,c("dim 2", "dim 3")])));
	colnames(cosDistance_dim2_vs_dim3) <- statsInputDf$labels_full;
	rownames(cosDistance_dim2_vs_dim3) <- statsInputDf$labels_full;

	cosDistance_dim1_vs_dim2 <- as.matrix(cosDistance_dim1_vs_dim2[varsNames_considered_rows,varsNames_catVarName]);
	colnames(cosDistance_dim1_vs_dim2) <- varsNames_catVarName
	cosDistance_dim1_vs_dim3 <- as.matrix(cosDistance_dim1_vs_dim3[varsNames_considered_rows,varsNames_catVarName]);
	colnames(cosDistance_dim1_vs_dim2) <- varsNames_catVarName
	cosDistance_dim2_vs_dim3 <- as.matrix(cosDistance_dim2_vs_dim3[varsNames_considered_rows,varsNames_catVarName]);
	colnames(cosDistance_dim1_vs_dim2) <- varsNames_catVarName
	cosDistanceTotal <- as.matrix(cosDistanceTotal[varsNames_considered_rows,varsNames_catVarName]);
	colnames(cosDistance_dim1_vs_dim2) <- varsNames_catVarName

	#save the full results
	resALL_CosDistances<- list(cosD_dim1_dim2=cosDistance_dim1_vs_dim2,cosD_dim1_dim3=cosDistance_dim1_vs_dim3,
		cosD_dim2_dim3=cosDistance_dim2_vs_dim3,cosD_All_dims=cosDistanceTotal);
	write.xlsx(resALL_CosDistances, file = paste0(wd, f_pheno_discrimination, f_results,f_cosineDistance,
		f_similarityThreshold,nameExperiment,"_",catVarName,"_cosSim_all_cosDistances_calculated.xlsx"),colNames = TRUE,rowNames = TRUE,);

	####################################################################
		#run the fucntion to extract the avriables with a % of cos simmilarity higher or equal to the threshold set
	####################################################################
	extract_cosSimilarity_vars_using_threshold.function(cosDistance_dim1_vs_dim2, catVarName, similarityThreshold,f_plot,"dim1_vs_dim2",nameExperiment);
	extract_cosSimilarity_vars_using_threshold.function(cosDistance_dim1_vs_dim3, catVarName, similarityThreshold,f_plot,"dim1_vs_dim3",nameExperiment);
	extract_cosSimilarity_vars_using_threshold.function(cosDistance_dim2_vs_dim3, catVarName, similarityThreshold,f_plot,"dim2_vs_dim3",nameExperiment);
	extract_cosSimilarity_vars_using_threshold.function(cosDistanceTotal, catVarName, similarityThreshold,f_plot,"all_dims",nameExperiment);
 };



#' Function for classification 
#'
#' @name glm_rf_classifiers_CatVarEffects.function 
#'
#' @param data input data can be 
#' For the full model, "phenoSocial_4analysis"
#' For the sel30 model, "parameters_more_30perc_informative_{nameVariable}"
#' For the sel1 model, "phenoAll_4analysis_woHighCor_sel1"
#' 
#' @param nameData name to assignate eEx:
#' For the full model, paste0(name,"_fullMdl")
#' For the sel30 model, paste0(name,"_sel30perc")
#' For the sel1 model, paste0(name,"_sel1Mdl")
#'
#' @param angle1 Default 90.The angles parameters represent the angles to 
#' write the variable names labels inside the plots,adjust if needed,
#' angle1 is for the glm plot

#' @param angle2 Default 90.The angles parameters represent the angles to 
#' write the variable names labels inside the plots, angle2 is for the RF plot.
#'
#' @param catVaMFA name of the independent variable used to run the mfa.function
#'
#' @param glmNet_grid keep it to normal, change to tuned when an error arise saying too many NAs were computed
#'
#' @param catVarName name of the categorical variable you want to test the contribution here.
#' 
#' @param nameModel Short naming for the model you are using.
#' Ex: "fullMdl" for the model including all variables;"Sel1" for the model not including the
#' highly correlated variables, "Sel30" for the model including only the variables contributing
#' more thana 30%
#' 
#' @return  calculation of cllassifiers.
#'
#' @examples
#'```
#' # Sel30 model (model most discriminative variables, 
#' contributing more than 30% in PCA)
#' nameModel <- "Sel30"
#' data.tmp <- parameters_more_30perc_informative_Genotype
#' glm_rf_classifiers_CatVarEffects.function(data.tmp,paste0(name,"_sel30perc"),90,90,"Genotype","Genotype","normal",nameModel);
#' Full model: all variables
#' nameModel <- "fullMdl"
#' data.tmp <- phenoSocial_4analysis
#' glm_rf_classifiers_CatVarEffects.function(data.tmp,paste0(name,"_fullMdl"),90,90,"Genotype","Genotype","normal",nameModel);
#' # Sel1 model (model without highly correlated variables)
#' nameModel <- "Sel1"
#' data.tmp <- phenoAll_4analysis_woHighCor_sel1
#' glm_rf_classifiers_CatVarEffects.function(data.tmp,paste0(name,"_sel1Mdl"),90,90,"Genotype","Genotype","normal",nameModel);
#'```
#' @export
glm_rf_classifiers_CatVarEffects.function <- function(data,nameData,angle1,angle2,catVaMFA,catVarName,glmNet_grid,nameModel){
	# description: classifiers algorithms to identify the main discriminative 
	# explanatory variables some specially selected to identify variables 
	# affected by gene dosage as GLM, multiGLM and GLMNet. 
	# And by an supervised method Random Forest or the unsupervised version of RF
	#NOTE:
	# The full_model is the model with all data : input file: phenoAll_4analysis_woHighCor
	# The model_selected_vars is the full model without the highly correlated variables

	# NOTE2 : to run the function glm_rf_classifiers_CatVarEffects.function you need to set several parameters
	#######################################################################################################
	## Function form: glm_rf_classifiers_CatVarEffects.function(data,nameData,angle1,angle2,catVaMFA,catVarName,glmNet_grid)
	### Parameters of the function between brackets: glm_rf_classifiers.function(data,nameData,angle1,angle2,glmNet_grid)
	### 1) "data" is the input dataframe containing all the variables info they begin by "phenoAll_4analysis_*"
	### 2) "nameData" is the name you want to give to the outputs of this function, as I will run the function 
	### with to input dataframes the one corresponding to all the model variables and another corresponding 
	### to the selected model variables I decided to use "full_model" and "sel_model".
	### 3) The angles parameters represent the angles to write the variable names labels inside the plots,
	###    adjust if needed angle1 is for the glm plot, angle2 is for the RF plot.
	### 4) glmNet_grid: keep it to normal, change to tuned when an error arise saying too many NAs were computed


    load( file=paste0(wd, f_pheno_discrimination,f_results, f_Rdata,name,"_",nameModel,"_",catVaMFA,"_Robject_for_plots.RData")) #namesVariablesDf, df1, df2, df3,df4, df5, df6, MFA1, MFA2, MFA3,nbVarF, nbVarTotalF,
    options(warn=-1)
	# Creating the folders for storing the data
	#######################################
	f_var_importance <- 'importance_variables/';
	f_varCatName <- paste0(gsub("^.* ","",catVarName), "/");
	catVarNameSaving <- gsub("^.* ","",catVarName);
	dir.create(file.path(wd, f_pheno_discrimination,f_results,f_var_importance), showWarnings = F);
	dir.create(file.path(wd, f_pheno_discrimination,f_results,f_var_importance,f_varCatName), showWarnings = F);
	###################################################
	
	data.tmp <- data;
	num_factors_outcome <- length(unique(data.tmp[,catVarName]));
	print(paste0("number of factors :", num_factors_outcome));
	#function to avoid having variance equal 0, at least 3 vals diffferent in all the data
	length_unique <- apply(data.tmp, 2, function(x) length(unique(x)))[-c(1:2)];
	which(length_unique[length_unique <2]=="TRUE");
    drops <- names(which((length_unique <3)=="TRUE"));
    data.tmp<- data.tmp[ , !(names(data.tmp) %in% drops)];
    varCatPosition <- which((colnames(data.tmp)==catVarName)==TRUE);
	data.tmp[,varCatPosition] <- factor(data.tmp[,varCatPosition], levels=unique(data.tmp[,varCatPosition]));

	col_idx <-colnames(data.tmp)[varCatPosition];
	data.tmp<- data.tmp %>% select(col_idx, everything());
	names_data.tmp <- colnames(data.tmp);
	colnames(data.tmp) <- paste0("var_", 1:ncol(data.tmp));

	if (num_factors_outcome==2) {

		####################################################
		# A1 step: Run the glm classifier
		#        <Generalized linear model > 
		####################################################
		
		# Creating the control
		#control <- trainControl(method="repeatedcv", number=1000, repeats=100);
		control <- trainControl(method = "cv", number = 10, # k-folds CV with k=10
                              classProbs = TRUE,
                              savePredictions = TRUE,
                              summaryFunction = multiClassSummary);# save predictions for ROC

	    #print(control_glm)
		# Running the modeldata=data.tmp,
		####data.data.tmp.cols <-data.tmp[,2:ncol(data.tmp)]
		####modelGlm <- caret::train(x=data.df.cols, y=data.tmp$var_1,  method="glm", preProcess="scale", trControl=control); #Generalized Linear Model 
		modelGlm <- caret::train(var_1~., data=data.tmp, method="glm", preProcess="scale", trControl=control); #Generalized Linear Model 
		#length( modelglm$finalModel$coefficients)-1 #nb selected variables
		assign(paste0(nameData,"_",catVarName,"_modelGlm"),modelGlm,.GlobalEnv);
	    save(modelGlm, file=paste0(wd,f_pheno_discrimination,f_results,f_var_importance,f_varCatName, nameData,"_modelGlm.RData"));
	    
	} else {
		
		####################################################
		# A2.A step: Run the multinom classifier 
		# < multinomial log-linear models via neural networks. > 
		####################################################
		# NOTE:
		# See the cran package for more info about the function is a multinomial glm
	
		model_multiGlm=multinom(var_1~., data=data.tmp); 

		####################################################
		# A2B step: Run the glmnet classifier
		#        < Elastic net algorithm > 
		####################################################
		# NOTE:
		# If the outcome variable we are predicting has more than two factors, we cannot use glm
		# to run a logistic regression. Instead, we will use glmnet to run an elastic net algorithm
		#  that can handle more factors in our outcome variable.
		# https://neurospection.netlify.app/post/machine-learning-basics-with-caret/
		# Creating the control
		control <- trainControl(method = "cv", number = 100, # k-folds CV with k=10
                              classProbs = TRUE,
                              savePredictions = TRUE,
                              summaryFunction = multiClassSummary);# save predictions for ROC

	#print(control_glm)
		if (glmNet_grid=="normal") {
	
		modelGlmNet <- train(var_1 ~. , 
	                      data = data.tmp, 
	                      method = "glmnet",
	                      trControl = control, 
	                      #tuneGrid = tune.grid,
	                      preProcess = c("center", "scale"),
	                      preProc = c("center"),
	                      standardize=FALSE);
		};
		if (glmNet_grid=="tuned") {
			
			tune.grid<- expand.grid(
				alpha = 0:1,
				lambda = seq(0.00001, 1, length = 1000));
	    
			modelGlmNet <- train(var_1 ~. , 
			data = data.tmp, 
			method = "glmnet",
			trControl = control, 
			tuneGrid = tune.grid,
			preProcess = c("center", "scale"),
			preProc = c("center"),
			standardize=FALSE);

			};

		save(modelGlmNet, file=paste0(wd,f_pheno_discrimination,f_results,f_var_importance,f_varCatName, nameData,"_modelGlmNet.RData"));
		#load(file=paste0(wd,f_pheno_discrimination,f_results,f_var_importance,f_varCatName, nameData,"_modelGlmNet.RData"));
		####################@ uncomment those lines to see further stats  ######################
		########################################################################################
		#glmnet.info <- getModelInfo("glmnet")
		#glmnet.info$glmnet$parameters
		#tune.grid <- expand.grid(alpha = modelGlmNet$finalModel$tuneValue$alpha,
	    #                     lambda = modelGlmNet$finalModel$tuneValue$lambda)
	    #
		#glmnet.predict <- predict(modelGlmNet, data.tmp) # Predict values in the testing set
		#postResample(glmnet.predict, data.tmp$var_1) # the accuracy of the model
		#confusionMatrix(glmnet.predict, data.tmp$var_1) 


		# Compute R^2 from true and predicted values 
		#The above output shows that the RMSE and R-squared values for the ridge regression model on the training data
		#eval_results <- function(true, predicted, df) {
		#	predicted <- as.numeric(gsub("homo","0", gsub("het","1", gsub("wt","2",predicted))));
		#	true <- as.numeric(gsub("homo","0", gsub("het","1", gsub("wt","2",true))));
		#  	SSE <- sum((predicted - true)^2)
		# 	SST <- sum((true - mean(true))^2)
		# 	R_square <- 1 - SSE / SST
		# 	RMSE = sqrt(SSE/nrow(df) )
		#	#glmnet_coefs(modelGlmNet, data.tmp[,-1], data.tmp$var_1, s = modelGlmNet$finalModel$lambdaOpt) #, prepend = "", ...
        #
	    #
		#	# Model performance metrics
		#	data.frame(RMSE = RMSE, Rsquare = R_square)
		#};
		#
		#res_Eval<-eval_results(data.tmp$var_1, glmnet.predict, data.tmp) 
		#glmnet_coefs(modelGlmNet, data.tmp[,-1], data.tmp$var_1, s = modelGlmNet$finalModel$lambdaOpt) #, prepend = "", ...
		########################################################################################

		assign(paste0(nameData,"_",catVarNameSaving,"_modelGlmNet"),modelGlmNet,.GlobalEnv);
		assign(paste0(nameData,"_",catVarNameSaving,"_model_multiGlm"),model_multiGlm,.GlobalEnv);


	};
	####################################################
	# B step: Run the RF classifier
	#            <Random forest > 
	####################################################

	# Creating the control

	#control_rf <- rfeControl(functions=rfFuncs, method="cv", number=1000, repeats=100);
	control <- trainControl(method = "cv", number = 10, # k-folds CV with k=10
                              classProbs = TRUE,
                              savePredictions = TRUE,
                              summaryFunction = multiClassSummary)# save predictions for ROC

	# Running the model
	modelRf <- caret::train(var_1~., data=data.tmp, method="rf", preProcess="scale", trControl=control,importance = TRUE); #RF
	####Conclusion: we go down in accuracy if we use sex:genotype all combinations 
	####modelRf2 <- caret::train(var_26~., data=data.tmp[,c(3:nbcols)], method="rf", preProcess="scale", trControl=control,importance = TRUE); #RF
	#assign(paste0(nameData,"_modelRf"),modelRf,.GlobalEnv);
    save(modelRf, file=paste0(wd,f_pheno_discrimination,f_results,f_var_importance,f_varCatName,nameData,"_modelRf.RData"));

	#######################################
	# Identify  the variables contributing 
	# more to the genotype discrimination 
	# and their scaled contribution
	#######################################
	scaled_importance_variables.function <- function(dataInput_4Function,classifier,nameClasifier,names_categories,num_factors_outcome){
		#README: classifierName can be ""glm"/ multiglm for more than 2 genotypes, glmNet for multiple co,parisons
		#  between each genotype more than 2 or "rf".
		#https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r
		#https://www.datacamp.com/community/tutorials/tutorial-ridge-lasso-elastic-net
		#https://statisticaloddsandends.wordpress.com/2018/11/15/a-deep-dive-into-glmnet-standardize/
		#http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/153-penalized-regression-essentials-ridge-lasso-elastic-net/
		#http://rstudio-pubs-static.s3.amazonaws.com/251240_12a8ecea8e144fada41120ddcf52b116.html#introducing-glmnet
		if (nameClasifier=="glm"){
	
			names_categories[1] <- paste0(catVarNameSaving);
			classifier_names <- data.frame(VariableName=names_categories,Variable=paste0("var_",1:(length(names_categories))));
			classifier$coefnames <-classifier_names$VariableName;
			classifiervarImp <- as.data.frame(varImp(classifier)$importance);
			classifiervarImp$Variable <- sub("^(\\D*\\d+).*", "\\1", rownames(classifiervarImp));
			classifiervarImp$classifier_2ndName <- gsub("var_.","",gsub("var_..$","",gsub("var_.$","",rownames(classifiervarImp))));
			classifiervarImp <- left_join(classifiervarImp, classifier_names,by="Variable") ;
			classifiervarImp$Variable <- gsub("_$","",paste0(classifiervarImp$VariableName, "_", classifiervarImp$classifier_2ndName));
		    drop <- c("classifier_2ndName","VariableName");
		    classifiervarImp<- classifiervarImp[ , !(colnames(classifiervarImp) %in% drop)]
            classifiervarImp$Variable<- factor(classifiervarImp$Variable,levels=unique(classifiervarImp$Variable));
            assign(paste0(nameClasifier,"_varImp"),classifiervarImp, envir=parent.frame());
            assign(nameClasifier,classifier, envir=parent.frame());
		
		} else if (nameClasifier=="glmNet") {

			classifier_name <- names_categories[-1];	#	
			classifier_names <- data.frame(VariableName=names_categories,Variable=paste0("var_",1:(length(names_categories))));
			#classifier_names <- data.frame(VariableName=c(classifier_name[1],classifier_name),Variable=paste0("var_",1:(length(classifier_name)+1)));
			classifiervarImp <- as.data.frame(varImp(classifier)$importance);
			classifiervarImp$Variable <- sub("^(\\D*\\d+).*", "\\1", rownames(classifiervarImp));
			classifiervarImp$classifier_2ndName <- gsub("var_.","",gsub("var_..$","",gsub("var_.$","",rownames(classifiervarImp))));
			classifiervarImp <- left_join(classifiervarImp, classifier_names,by="Variable") ;
			classifiervarImp$Variable <- gsub("_$","",paste0(classifiervarImp$VariableName, "_", classifiervarImp$classifier_2ndName));
		    drop <- c("classifier_2ndName","VariableName");
		    classifiervarImp<- classifiervarImp[ , !(colnames(classifiervarImp) %in% drop)]
            classifiervarImp$Variable<- factor(classifiervarImp$Variable,levels=unique(classifiervarImp$Variable));
            assign(paste0(nameClasifier,"_varImp"),classifiervarImp, envir=parent.frame());
            assign(nameClasifier,classifier, envir=parent.frame());

        } else if (nameClasifier=="rf") {
		
			classifier_name <- names_categories[-1];	#	
			classifier_names <- data.frame(VariableName=names_categories,Variable=paste0("var_",1:(length(names_categories))));
			#classifier_names <- data.frame(VariableName=c(classifier_name[1],classifier_name),Variable=paste0("var_",1:(length(classifier_name)+1)));
			classifiervarImp <- as.data.frame(varImp(classifier)$importance);
			classifiervarImp$Variable <- sub("^(\\D*\\d+).*", "\\1", rownames(classifiervarImp));
			classifiervarImp$classifier_2ndName <- gsub("var_.","",gsub("var_..$","",gsub("var_.$","",rownames(classifiervarImp))));
			classifiervarImp <- left_join(classifiervarImp, classifier_names,by="Variable") ;
			classifiervarImp$Variable <- gsub("_$","",paste0(classifiervarImp$VariableName, "_", classifiervarImp$classifier_2ndName));
		    drop <- c("classifier_2ndName","VariableName");
		    classifiervarImp<- classifiervarImp[ , !(colnames(classifiervarImp) %in% drop)]
            classifiervarImp$Variable<- factor(classifiervarImp$Variable,levels=unique(classifiervarImp$Variable));

            assign(paste0(nameClasifier,"_varImp"),classifiervarImp, envir=parent.frame());
            assign(nameClasifier,classifier, envir=parent.frame());

        } else if (nameClasifier=="multiGlm") {

			classifier_name <- names_categories[-1];	#	
			classifier_names <- data.frame(VariableName=names_categories,Variable=paste0("var_",1:(length(names_categories))));
			#classifier_names <- data.frame(VariableName=c(classifier_name[1],classifier_name),Variable=paste0("var_",1:(length(classifier_name)+1)));
			classifiervarImp <- as.data.frame(varImp(classifier));
			classifiervarImp$Variable <- sub("^(\\D*\\d+).*", "\\1", rownames(classifiervarImp));
			classifiervarImp$classifier_2ndName <- gsub("var_.","",gsub("var_..$","",gsub("var_.$","",rownames(classifiervarImp))));
			classifiervarImp <- left_join(classifiervarImp, classifier_names,by="Variable") ;
			classifiervarImp$Variable <- gsub("_$","",paste0(classifiervarImp$VariableName, "_", classifiervarImp$classifier_2ndName));
		    drop <- c("classifier_2ndName","VariableName");
		    classifiervarImp<- classifiervarImp[ , !(colnames(classifiervarImp) %in% drop)]
            classifiervarImp$Variable<- factor(classifiervarImp$Variable,levels=unique(classifiervarImp$Variable));
            assign(paste0(nameClasifier,"_varImp"),classifiervarImp, envir=parent.frame());
            assign(nameClasifier,classifier, envir=parent.frame());


		} else {
			print("Error, the statistical method asked for has not being defined, check your labelling :)")

		};   
	};
	 
 	if (num_factors_outcome==2) {
   
    	print("factors <=  2 beginning the glm and rf analyses")
	 	# run the annidated function
		scaled_importance_variables.function(data.tmp,modelGlm, "glm",names_data.tmp,num_factors_outcome);
		print("glm: ok ;)");
		scaled_importance_variables.function(data.tmp,modelRf, "rf",names_data.tmp,num_factors_outcome);
		print("rf: ok ;)");
	    print("1) The RF/Glm models were succesfully created :) ");

	    #save the vars into the global env to be able to access them
		assign(paste0(nameData,"_", catVarNameSaving,"_rf_varImp"),rf_varImp,.GlobalEnv);	
		assign(paste0(nameData,"_", catVarNameSaving,"_glm_varImp"),glm_varImp,.GlobalEnv);	
		assign(paste0(nameData,"_", catVarNameSaving,"_rf_classifier"),modelRf,.GlobalEnv);
		assign(paste0(nameData,"_", catVarNameSaving,"_glm_classifier"),modelGlm,.GlobalEnv);

		#store the results in a list. Each element of the list will be copied into a sheet in xlsx
		classifiers_importance <- list(glm=glm_varImp,rf=rf_varImp);
		openxlsx::write.xlsx(classifiers_importance, file = paste0(wd,f_pheno_discrimination,f_results,f_var_importance,f_varCatName,nameData, "_scaled_importance_by_glm_and_rf.xlsx"));

		#######################################
		### Visualization classifier glm & rf
		#######################################
	    options(warn=0);
		
		#for GLM
		glm_varImp$Variable <- factor(glm_varImp$Variable,levels=unique(glm_varImp$Variable))
		plot_glm <- ggplot(glm_varImp) + geom_point(aes(x=Variable, y=Overall,size=5),shape=16); 
		plot_glm <- plot_glm + scale_y_continuous("scaled importance",breaks=c(0,25,50,75,100),limits=c(0,120)) + 
		theme(axis.text.x = element_text(angle = angle1, hjust = 1, vjust = 1,size=7),axis.title=element_text(size=13,face="bold"),legend.position="none",plot.title = element_text(size=6))+
		scale_x_discrete("",expand = waiver()) +labs( title=paste("Vars importance: ",nameData,", ", nbVarF," vars: ",catVarNameSaving," : glm" ), caption = "Source: Y.Herault team @IGBMC") +   
		guides(size=FALSE)+ coord_flip() + theme(axis.text.x = element_text(size = 15,angle=0),axis.text.y=element_text(size=12),plot.margin = unit(c(1,1,1,1), "cm")); 
		
		#for RF
		if (ncol(rf_varImp)==2){
			colnames(rf_varImp)[1] <- "Overall"
			rf_varImp$Variable <- factor(rf_varImp$Variable,levels=unique(rf_varImp$Variable))
			plot_rf<- ggplot(rf_varImp) + geom_point(aes(x=Variable, y=Overall,size=5,col="Black"),shape=16); 
			plot_rf<- plot_rf+ scale_y_continuous("scaled importance",breaks=c(0,25,50,75,100),limits=c(0,120)) + 
			theme(axis.text.x = element_text(angle = angle2, hjust = 1, vjust = 1,size=7),axis.title=element_text(size=13,face="bold"),legend.position="none",plot.title = element_text(size=6))+
			scale_x_discrete("") +labs( title=paste("Vars importance: ",nameData, ", ", nbVarF," vars:",catVarNameSaving," : RF" ), caption = "Source: Y.Herault team @IGBMC") +   
			guides(size=FALSE)+ coord_flip() + theme(axis.text.x = element_text(size = 15,angle=0),axis.text.y=element_text(size=12),plot.margin = unit(c(1,1,1,1), "cm"))+
			scale_color_manual(values="Black"); 
		} else{
			
			if(rf_varImp[1,1]==rf_varImp[1,2] && rf_varImp[2,1]==rf_varImp[2,2] ){
				rf_varImp_varImp_input_4_multiPlot <- rf_varImp[,-1];
				colnames(rf_varImp_varImp_input_4_multiPlot)[1] <-"Overall";
				rf_varImp_varImp_input_4_multiPlot$Variable <- factor(rf_varImp_varImp_input_4_multiPlot$Variable,levels=unique(rf_varImp_varImp_input_4_multiPlot$Variable))
				plot_rf<- ggplot(rf_varImp_varImp_input_4_multiPlot) + geom_point(aes(x=Variable, y=Overall,size=5,col="Black"),shape=16); 
				plot_rf<- plot_rf+ scale_y_continuous("scaled importance",breaks=c(0,25,50,75,100),limits=c(0,120)) + 
				theme(axis.text.x = element_text(angle = angle2, hjust = 1, vjust = 1,size=7),axis.title=element_text(size=13,face="bold"),legend.position="none",plot.title = element_text(size=6))+
				scale_x_discrete("") +labs( title=paste("Vars importance: ",nameData, ", ", nbVarF," vars:",catVarNameSaving," : RF" ), caption = "Source: Y.Herault team @IGBMC") +   
				guides(size=FALSE)+ coord_flip() + theme(axis.text.x = element_text(size = 15,angle=0),axis.text.y=element_text(size=12),plot.margin = unit(c(1,1,1,1), "cm"))+	
				scale_color_manual(values="Black"); 
; 
			} else {
				rf_varImp_varImp_input_4_multiPlot <- melt(rf_varImp, id.vars="Variable");
				rf_varImp_varImp_input_4_multiPlot<- left_join(rf_varImp_varImp_input_4_multiPlot,namesVariablesDf,by="Variable");

				colnames(rf_varImp_varImp_input_4_multiPlot)[which(colnames(rf_varImp_varImp_input_4_multiPlot)=="variable")] <- "Group";
				colnames(rf_varImp_varImp_input_4_multiPlot)[which(colnames(rf_varImp_varImp_input_4_multiPlot)=="value")] <- "Overall";
				rf_varImp_varImp_input_4_multiPlot$Group <- gsub("_.\\.","_",rf_varImp_varImp_input_4_multiPlot$Group);


				namesVariablesDfmod <- namesVariablesDf[,-1];
				colnames(namesVariablesDfmod) <- c("colPaletteGroup","Group")
				rf_varImp_varImp_input_4_multiPlot<- left_join(rf_varImp_varImp_input_4_multiPlot,namesVariablesDfmod,by="Group");
				rf_varImp_varImp_input_4_multiPlot$Variable <- factor(rf_varImp_varImp_input_4_multiPlot$Variable,levels=unique(rf_varImp_varImp_input_4_multiPlot$Variable))
				rf_varImp_varImp_input_4_multiPlot$Group <- factor(rf_varImp_varImp_input_4_multiPlot$Group,levels=unique(rf_varImp_varImp_input_4_multiPlot$Group))


				#rf_varImp_varImp_input_4_multiPlot$colPaletteGroup<- namesVariablesDf[which(namesVariablesDf$Variable  %in%  rf_varImp_varImp_input_4_multiPlot[,'Group']), 'colPalette']
				plot_rf <- ggplot(rf_varImp_varImp_input_4_multiPlot) + geom_point(aes(x=Variable, y=Overall,col=Group,size=5),shape=16)+ 
				scale_color_manual(values=rf_varImp_varImp_input_4_multiPlot$colPaletteGroup,labels=rf_varImp_varImp_input_4_multiPlot$Group); 
				plot_rf <- plot_rf + scale_y_continuous("scaled importance",breaks=c(0,25,50,75,100),limits=c(0,120)) + 
				theme(axis.text.x = element_text(angle = angle1, hjust = 1, vjust = 1,size=7),axis.title=element_text(size=13,face="bold"),plot.title = element_text(size=6))+
				scale_x_discrete("",expand = waiver()) +labs( title=paste("Vars importance: ",nameData,", ", nbVarF," vars: ",catVarNameSaving," : RF" ), caption = "Source: Y.Herault team @IGBMC") +   
				guides(size=FALSE)+ coord_flip() + theme(axis.text.x = element_text(size = 15,angle=0),axis.text.y=element_text(size=12),plot.margin = unit(c(1,1,1,1), "cm"),,legend.position="bottom"); 		
			};
		};

		assign(paste0("plot_glm_",catVarNameSaving,"_",nameData),plot_glm,.GlobalEnv);
		assign(paste0("plot_rf_",catVarNameSaving,"_",nameData),plot_rf,.GlobalEnv);    
   		print("All jobs done, pass to the next cell :)");
    	
    } else {
    	print("factors > than 2 beginning the analysis")
   	
		scaled_importance_variables.function(data.tmp,model_multiGlm, "multiGlm",names_data.tmp,num_factors_outcome);
		print("multiGlm: ok ;)")
		scaled_importance_variables.function(data.tmp,modelGlmNet, "glmNet",names_data.tmp,num_factors_outcome);
		print("glmNet: ok ;)")
		scaled_importance_variables.function(data.tmp,modelRf, "rf",names_data.tmp,num_factors_outcome);
		print("rf: ok ;)")
	    print("1) The RF/multiGlm/glmNet models were succesfully created :) ");

	    #save the vars into the global env to be able to access them 
		assign(paste0(nameData,"_",catVarNameSaving,"_multiGlm_varImp"),multiGlm_varImp,.GlobalEnv);	
		assign(paste0(nameData,"_",catVarNameSaving,"_multiGlm_classifier"),model_multiGlm ,.GlobalEnv);
		assign(paste0(nameData,"_",catVarNameSaving,"_glmNet_varImp"),glmNet_varImp,.GlobalEnv);	
		assign(paste0(nameData,"_",catVarNameSaving,"_glmNet_classifier"),modelGlmNet,.GlobalEnv);
		assign(paste0(nameData,"_",catVarNameSaving,"_rf_varImp"),rf_varImp,.GlobalEnv);	
		assign(paste0(nameData,"_",catVarNameSaving,"_rf_classifier"),modelRf,.GlobalEnv);

		#store the results in a list. Each element of the list will be copied into a sheet in xlsx
		classifiers_importance <- list(glmNet=glmNet_varImp, multiGlm=multiGlm_varImp,rf=rf_varImp);
		openxlsx::write.xlsx(classifiers_importance, file = paste0(wd,f_pheno_discrimination,f_results,f_var_importance,f_varCatName,nameData, "_scaled_importance_by_multiGlm_glmNet_and_rf.xlsx"));
		
		#######################################
		### Visualization classifier glm & rf
		#######################################

		## Preprocesing step: formatting the VarImp dd to plot for each factor value of the categorical variables 
		##  the most important variables in diffeent colours 

		rf_varImp_varImp_input_4_multiPlot <- melt(rf_varImp, id.vars="Variable");
		colnames(rf_varImp_varImp_input_4_multiPlot)[which(colnames(rf_varImp_varImp_input_4_multiPlot)=="variable")] <- "Group";
		colnames(rf_varImp_varImp_input_4_multiPlot)[which(colnames(rf_varImp_varImp_input_4_multiPlot)=="value")] <- "Overall";
		rf_varImp_varImp_input_4_multiPlot<- left_join(rf_varImp_varImp_input_4_multiPlot,namesVariablesDf,by="Variable");
		rf_varImp_varImp_input_4_multiPlot$Group <- gsub("_.\\.","_",rf_varImp_varImp_input_4_multiPlot$Group);
		namesVariablesDfmod <- namesVariablesDf[,-1];
		colnames(namesVariablesDfmod) <- c("colPaletteGroup","Group")
		print(namesVariablesDfmod)
		rf_varImp_varImp_input_4_multiPlot<- left_join(rf_varImp_varImp_input_4_multiPlot,namesVariablesDfmod,by="Group");
		rf_varImp_varImp_input_4_multiPlot$Variable <- factor(rf_varImp_varImp_input_4_multiPlot$Variable,levels=unique(rf_varImp_varImp_input_4_multiPlot$Variable))
		rf_varImp_varImp_input_4_multiPlot$Group <- factor(rf_varImp_varImp_input_4_multiPlot$Group,levels=unique(rf_varImp_varImp_input_4_multiPlot$Group))
		
		glmNet_varImp_input_4_multiPlot <- melt(glmNet_varImp, id.vars="Variable");
		colnames(glmNet_varImp_input_4_multiPlot)[which(colnames(glmNet_varImp_input_4_multiPlot)=="variable")] <- "Group";
		colnames(glmNet_varImp_input_4_multiPlot)[which(colnames(glmNet_varImp_input_4_multiPlot)=="value")] <- "Overall";
		glmNet_varImp_input_4_multiPlot<- left_join(glmNet_varImp_input_4_multiPlot,namesVariablesDf,by="Variable");
		glmNet_varImp_input_4_multiPlot$Group <- gsub("_.\\.","_",glmNet_varImp_input_4_multiPlot$Group);
		glmNet_varImp_input_4_multiPlot<- left_join(glmNet_varImp_input_4_multiPlot,namesVariablesDfmod,by="Group");
		glmNet_varImp_input_4_multiPlot$Variable <- factor(glmNet_varImp_input_4_multiPlot$Variable,levels=unique(glmNet_varImp_input_4_multiPlot$Variable))
		glmNet_varImp_input_4_multiPlot$Group <- factor(glmNet_varImp_input_4_multiPlot$Group,levels=unique(glmNet_varImp_input_4_multiPlot$Group))



		#plot functions:
		maxVal <- as.numeric(max(glmNet_varImp_input_4_multiPlot$Overall)+10);
		plot_glmNet <- ggplot(glmNet_varImp_input_4_multiPlot) + geom_point(aes(x=Variable, y=Overall,col=Group,size=5),shape=16)+ 		
		scale_color_manual(values=unique(glmNet_varImp_input_4_multiPlot$colPaletteGroup),labels=unique(glmNet_varImp_input_4_multiPlot$Group));  
		plot_glmNet <- plot_glmNet + scale_y_continuous("scaled importance",breaks=c(0,25,50,75,100),limits=c(0,120)) + 
		theme(axis.text.x = element_text(angle = angle1, hjust = 1, vjust = 1,size=7),axis.title=element_text(size=13,face="bold"),plot.title = element_text(size=6))+
		scale_x_discrete("",expand = waiver()) +labs( title=paste("Vars importance: ",nameData,", ", nbVarF," vars: ",catVarNameSaving," : glmNet" ), caption = "Source: Y.Herault team @IGBMC") +   
		guides(size=FALSE)+ coord_flip() + theme(axis.text.x = element_text(size = 15,angle=0),axis.text.y=element_text(size=12),plot.margin = unit(c(1,1,1,1), "cm"),legend.position="bottom")

		plot_rf <- ggplot(rf_varImp_varImp_input_4_multiPlot) + geom_point(aes(x=Variable, y=Overall,col=Group,size=5),shape=16)+
	    scale_color_manual(values=unique(rf_varImp_varImp_input_4_multiPlot$colPaletteGroup),labels=unique(rf_varImp_varImp_input_4_multiPlot$Group)); 
		plot_rf <- plot_rf + scale_y_continuous("scaled importance",breaks=c(0,25,50,75,100),limits=c(0,120)) + 
		theme(axis.text.x = element_text(angle = angle1, hjust = 1, vjust = 1,size=7),axis.title=element_text(size=13,face="bold"),plot.title = element_text(size=6))+
		scale_x_discrete("",expand = waiver()) +labs( title=paste("Vars importance: ",nameData,", ", nbVarF," vars: ",catVarNameSaving," : RF" ), caption = "Source: Y.Herault team @IGBMC") +   
		guides(size=FALSE)+ coord_flip() + theme(axis.text.x = element_text(size = 15,angle=0),axis.text.y=element_text(size=12),plot.margin = unit(c(1,1,1,1), "cm"),legend.position="bottom"); 

		plot_multiglm <- ggplot(multiGlm_varImp) + geom_point(aes(x=Variable, y=multiGlm_varImp[,1],size=5),shape=16); 
		plot_multiglm <- plot_multiglm + scale_y_continuous("scaled importance") + 
		theme(axis.text.x = element_text(angle = angle1, hjust = 1, vjust = 1,size=7),axis.title=element_text(size=13,face="bold"),plot.title = element_text(size=6))+
		scale_x_discrete("",expand = waiver()) +labs( title=paste("Vars importance: ",nameData,", ", nbVarF, " vars: ",catVarNameSaving," : multiGLM" ), caption = "Source: Y.Herault team @IGBMC") +   
		guides(size=FALSE)+ coord_flip() + theme(axis.text.x = element_text(size = 15,angle=0),axis.text.y=element_text(size=12),plot.margin = unit(c(1,1,1,1), "cm"));
		assign(paste0("plot_glmNet_",catVarNameSaving,"_",nameData),plot_glmNet,.GlobalEnv);
		assign(paste0("plot_multiglm_",catVarNameSaving,"_",nameData),plot_multiglm,.GlobalEnv);
		assign(paste0("plot_rf_",catVarNameSaving,"_",nameData),plot_rf,.GlobalEnv);
	    options(warn=0);
	
	    # Save all the plots generated
	    #grDevices::pdf(paste0(wd,f_pheno_discrimination, f_results,f_var_importance,nameData, "_variables_importance.pdf"));
	    #plot_glmNet;
	    #plot_multiglm;
	    #plot_rf;
	    #dev.off();
	    #print(" 2) The plotting finished  well too :)");
    
   	print("All jobs done, pass to the next cell :)");
    };
};



#' Parallel plot
#'
#' @name stats_calculation_for_a_df_with_observations.function 
#'
#' @param inputFile formatted as shown in vignette from the files:
#' For the full model, "phenoSocial_4analysis"
#' For the sel30 model, "parameters_more_30perc_informative_{nameVariable}"
#' For the sel1 model, "phenoAll_4analysis_woHighCor_sel1"
#' 
#' @param nameOutput Sufix Name specific from this project to add to the name of the output of the function
#' 
#' @param nbCatVars Number of categorical variables contained in dataInputPCA. Ex use 2 if Genotype and Sex  are present.
#' 
#' @return  calculation of cllassifiers.
#'
#' @examples
#' stats_calculation_for_a_df_with_observations.function(inputDf_parallele_Ind,"Genotype",1)
#'
#' @export
stats_calculation_for_a_df_with_observations.function <- function(inputFile,outputName,nbCatVars){
	# description: calculation of mean, stand error, CI to produce the 
	# formatted input needed to compute the parallel plot

	errorDf <- inputFile[,];
	nbCatVarsIndex<- nbCatVars+1
	catVarsName <- colnames(inputFile)[nbCatVars];
	errorDf[,c(nbCatVarsIndex:length(errorDf))] <- NA;
	rownames(errorDf)<- NULL;
	errorDf<- unique(errorDf);
	ci_upper_Df<- errorDf;
	ci_lower_Df<- errorDf;
	ci_mean_Df <-errorDf;
	sd_Df <-errorDf;
	for (i in nbCatVarsIndex:ncol(inputFile)) {
		df <- inputFile[,c(1,i)];
		for (j in 1:length(unique(df[,1]))) {
			dfj<- df[which(df$condition==unique(df$condition)[j]),];
			errorDf[j,i] <-qt(0.975,df=length(dfj[,2])-1)*sd(dfj[,2])/sqrt(length(dfj[,2]));
			ci_upper_Df[j,i] <-CI(dfj[,2], ci = 0.95)[1];
			ci_lower_Df[j,i] <-CI(dfj[,2], ci = 0.95)[3];
			ci_mean_Df[j,i] <-CI(dfj[,2], ci = 0.95)[2];
			sd_Df[j,i] <-sd(dfj[,2]);

			j=j+1;
		}
		i=i+1;
	}
	assign(paste0("error_",outputName),errorDf,.GlobalEnv);
	assign(paste0("CIupper_",outputName),ci_upper_Df,.GlobalEnv);
	assign(paste0("CIlower_",outputName),ci_lower_Df,.GlobalEnv);
	assign(paste0("sd_",outputName),ci_mean_Df,.GlobalEnv);
	assign(paste0("mean_",outputName),sd_Df,.GlobalEnv);

	colnames(errorDf)[nbCatVarsIndex:ncol(errorDf)] <- paste0("error ", colnames(errorDf)[nbCatVarsIndex:ncol(errorDf)] );
	colnames(ci_upper_Df)[nbCatVarsIndex:ncol(ci_upper_Df)] <- paste0("ci_upper ", colnames(ci_upper_Df)[nbCatVarsIndex:ncol(ci_upper_Df)] );
	colnames(ci_lower_Df)[nbCatVarsIndex:ncol(ci_lower_Df)] <- paste0("ci_lower ", colnames(ci_lower_Df)[nbCatVarsIndex:ncol(ci_lower_Df)] );
	colnames(ci_mean_Df)[nbCatVarsIndex:ncol(ci_mean_Df)] <- paste0("mean ", colnames(ci_mean_Df)[nbCatVarsIndex:ncol(ci_mean_Df)] );
	colnames(sd_Df)[nbCatVarsIndex:ncol(sd_Df)] <- paste0("sd ", colnames(sd_Df)[nbCatVarsIndex:ncol(sd_Df)] );

	pPlot_stats_Df <- left_join(ci_mean_Df,sd_Df,by=unique(c("condition",catVarsName)));
	pPlot_stats_Df <- left_join(pPlot_stats_Df,errorDf,by=unique(c("condition",catVarsName)));
	pPlot_stats_Df <- left_join(pPlot_stats_Df,ci_upper_Df,by=unique(c("condition",catVarsName)));
	pPlot_stats_Df <- left_join(pPlot_stats_Df,ci_lower_Df,by=unique(c("condition",catVarsName)));
	
	assign(paste0("inputDf_parallel_Sumup_",outputName),pPlot_stats_Df,.GlobalEnv);

};


