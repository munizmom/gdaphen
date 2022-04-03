#' Preformatting- Annonymization adn correlation
#'
#' @name condition_definition.function 
#'
#' @param InputDf Is  the phenoInput file , orginal file after imputation if is needed, 
#' of the data provided by the experimentator.
#'
#' @param conditionNameType "perRowNumber"  when we want to annonymize or run Gdaphen in more than 
#'  1 model or one model with different ages/treatments ,when the name labels are long, etc. 
#'  Instead we will use conditionNameType=Ind when is only one model, with one age:treatment data
#'	By default recommended conditionNameType=perRowNumber
#'
#' @param catVariables_initial_colNb Number of qualitative variables on InputDf, perform a 
#' head(InputDf) to get this number.
#' 
#' @return Several objects, the input for the full model phenoSocial_4analysis,
#' the imput for the selected model without highly corr variables, noted as Sel1, 
#' phenoAll_4analysis_woHighCor_sel1. Both with and without Genotype.
#' The correlationMatrix, and other parameters for downstream functions.
#' pheno: Is the InputDf, assigned to anew name to avoid touching the input data.
#'
#' @examples
#'	condition_definition.function(phenoInput,"perRowNumber");
#'	condition_definition.function(phenoInput,"Ind"); 
#' condition_definition.function(phenoInput,"Ind",5);
#' @export
condition_definition.function<- function(inputDf,conditionNameType,catVariables_initial_colNb){
	#####################################################################################################
	#DESCRIPTION: Function to annonymize data if it comes from human samples or a desing double blinded wants to be performed.
	# Its also works to rename animals with same genotype, sex and Ind number. 
	# This function also compute the correlation between variables.
	#########################
	#NOte: we will use by default conditionNameType="perRowNumber"  when we run Gdaphen in more than 
	# 1 model or one model with different ages/treatmetns etc. 
	# Instead we will use conditionNameType="Ind" when is only one model, with one age/treatment data
	#running examples:
	# > condition_definition.function(phenoInput,"perRowNumber");
	# > condition_definition.function(phenoInput,"Ind"); 
	#catVariables_initial_colNb, nb of columns that are categorical after the initial Ind column
	#####################################################################################################

	if (conditionNameType=="Ind"){
		inputDf$condition <- paste0(gsub("male","M_",gsub("female","F_",inputDf$`Sex:: Sex`)),inputDf$Genotype,"_",inputDf$Ind);

		dupIDs <- inputDf$condition[duplicated(inputDf$condition)];
		#renaming rownames if Ind, sex,Genotype conbinations are similar in any case
		############################@

		if (length(unique(dupIDs))>0){
			for (i in 1:length(unique(dupIDs))) {
					dupRows <- inputDf[which(inputDf$condition==unique(dupIDs)[i]),];
					dupRows$condition <- paste0(dupRows$condition, "_nb",1:nrow(dupRows));
					nonDupDf<- inputDf[-which(inputDf$condition==unique(dupIDs)[i]),];
					inputDf <- bind_rows(nonDupDf,dupRows);
					i <- i+1;
			};
		};
	rownames(inputDf) <- inputDf$condition;
	} else {
		inputDf$condition <- paste0(gsub("male","M_",gsub("female","F_",inputDf$`Sex:: Sex`)),inputDf$Genotype,"_",1:nrow(inputDf));
		rownames(inputDf) <- inputDf$condition;
	};

	colnames(inputDf) <- gsub(" $","",gsub("        $","",colnames(inputDf) ));
    drops <- c("Ind","condition");
    phenoSocial_4analysis<- inputDf[ , !(names(inputDf) %in% drops)];
	nb_vars <- ncol(phenoSocial_4analysis);

	# Preprocessing step: 
	# Removing the categorical variables with less than 2 diff values or the quantitative variables with less than 3 different values, 
	############################@
	catVarsDf<- phenoSocial_4analysis[,c(1:catVariables_initial_colNb)];
	catVar_selected <- colnames(Filter(function(x)(length(unique(x))>1), catVarsDf))
	catVariables_4_analysis_nb <- length(catVar_selected) #nb of cols of categorical variables where there are 2 or more different factors possible  

	quantitativeVarsDf <- phenoSocial_4analysis[,c(catVariables_initial_colNb:length(colnames(phenoSocial_4analysis)))];
	quantitativeVars_selected <- colnames(Filter(function(x)(length(unique(x))>2), quantitativeVarsDf))
	allVars_selected_names <- c(catVar_selected, quantitativeVars_selected) 

	#final input file:
	phenoSocial_4analysis <- phenoSocial_4analysis[,allVars_selected_names]

	####################################################
	# Step 1: remove variables highly correlated
	####################################################
	# step 1a:  calculate correlation matrix to remove variables highly correlated from the analysis 
	#parameter to set up:
	nb_varsF <- ncol(phenoSocial_4analysis);
	correlationMatrix <- cor(phenoSocial_4analysis[,c(as.numeric(catVariables_4_analysis_nb+1):nb_varsF)]); # wo qualitative variables

	# Step 1b: find attributes that are highly corrected (ideally >0.75);
	highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75);
	highlyCorrelated.Colnames <- colnames(correlationMatrix[,highlyCorrelated]);

	if (length(highlyCorrelated.Colnames) > 0) {
		HighCorrelatedVarsMatrix <- as.data.frame(correlationMatrix[highlyCorrelated.Colnames,highlyCorrelated.Colnames]);	  #saving the high corr variables in an excel
		HighCorrelatedVarsMatrix$variable <- rownames(HighCorrelatedVarsMatrix);
		HighCorrelatedVarsMatrixAll <- as.data.frame(correlationMatrix[,highlyCorrelated.Colnames]);
		HighCorrelatedVarsMatrixAll$variable <- rownames(HighCorrelatedVarsMatrixAll);
		corrRes <- list(CorrelationsCalculated=correlationMatrix,corVaribales075=highlyCorrelated.Colnames);

		highlyCorrelated.tmp <- correlationMatrix[,highlyCorrelated.Colnames]; 		#checking the identity of the correlated variables with each of the ones identified
		correlatedVars_listRes<- list();
		for (i in 1:ncol(highlyCorrelated.tmp)) {
			variableTested <- colnames(highlyCorrelated.tmp)[i];
			variableAllCorr <-highlyCorrelated.tmp[,i];
			variableAllCorrT <- abs(variableAllCorr)>=0.75;
			variablesCorrelated_with_first <- names(variableAllCorrT[variableAllCorrT == "TRUE" ]);
			variablesCorrelated_with_first_values <- as.data.frame(correlationMatrix[variablesCorrelated_with_first,variableTested]);
			colnames(variablesCorrelated_with_first_values) <- variableTested;
			correlatedDf<- as.data.frame(t(variablesCorrelated_with_first_values));
			colnames(correlatedDf) <- colnames(correlatedDf);
			correlatedDf$variableTested_for_correlation <- variableTested;
			correlatedVars_listRes[[i]]<- correlatedDf
			i <- i+1;
	};

	correlatedVars <- append(corrRes, correlatedVars_listRes);
	write.xlsx(correlatedVars, file = paste0(wd,f_pheno_discrimination, f_results, f_correlatedVars,name ,"_highlycorrelated_vars.xlsx"),colnames=T,rownames=T);
	assign("correlatedVars",correlatedVars ,.GlobalEnv);
	assign("highlyCorrelated",highlyCorrelated ,.GlobalEnv);
	assign("highlyCorrelated.Colnames",highlyCorrelated.Colnames ,.GlobalEnv);
	phenoAll_4analysis_woHighCor_sel1 <-phenoSocial_4analysis[, !(colnames(phenoSocial_4analysis) %in% highlyCorrelated.Colnames)]; #model4
	} else {
		phenoAll_4analysis_woHighCor_sel1 <-phenoSocial_4analysis;
	};

    dropsGeno <-"Genotype"
    phenoAll_4analysis_woHighCor_sel1woGenotype<- phenoAll_4analysis_woHighCor_sel1[ , !(names(phenoAll_4analysis_woHighCor_sel1) %in% dropsGeno)];
    phenoSocial_4analysis_selwoGenotype<- phenoSocial_4analysis[ , !(names(phenoSocial_4analysis) %in% dropsGeno)];

	assign("pheno",inputDf ,.GlobalEnv);
	assign("phenoSocial_4analysis",phenoSocial_4analysis ,.GlobalEnv);
	assign("nb_vars",nb_vars ,.GlobalEnv);
	assign("catVariables_initial_colNb",catVariables_initial_colNb ,.GlobalEnv);
	assign("allVars_selected_names",allVars_selected_names ,.GlobalEnv);
	assign("correlationMatrix",correlationMatrix ,.GlobalEnv);
	assign("phenoSocial_4analysis_selwoGenotype",phenoSocial_4analysis_selwoGenotype ,.GlobalEnv);
	assign("phenoAll_4analysis_woHighCor_sel1",phenoAll_4analysis_woHighCor_sel1 ,.GlobalEnv);
	assign("phenoAll_4analysis_woHighCor_sel1woGenotype",phenoAll_4analysis_woHighCor_sel1woGenotype ,.GlobalEnv);
	assign("nb_varsF",nb_varsF ,.GlobalEnv);
	assign("catVar_selected",catVar_selected ,.GlobalEnv);
};



#' Preformatting colour per sample.
#'
#' @name indSexCond.colour_definition.function 
#'
#'
#' @param pheno Is  the InputDf file , orginal file after imputation if is needed, 
#' of the data provided by the experimentator.
#'
#' @param catVar_selected name of the categorical variables that are selected for the analyses
#'
#' @param name_catVar_selected Rename if you want the columns names representing the categorical variable selected.
#' 
#' @return Dataframe with the name fo the samples, and the color asignation for each Individual point depending
#' on the categorical variables levels. Ex: depending on genotypes, sex, 
#' treatments. This function is advised to be modified by the user as you may have adifferent experimental 
#  design than us. The result of this function is formated along the vignette please take a look
#'
#' @examples
#' indSexCond.colour_definition.function(pheno, catVar_selected,catVar_selected);
#'
#' @export
indSexCond.colour_definition.function <- function(pheno,catVar_selected,name_catVar_selected){
	# description Function to annonymize data if it comes from human samples or a desing double blinded wants to be performed.
	# Its also works to rename animals with same genotype, sex and Ind number. 
	# This function also compute the correlation between variables.
	
	if (length(catVar_selected)>1) {
		indSexCond <- pheno[,catVar_selected];
		indSexCond$labels <- rownames(indSexCond);

	} else {
		indSexCond <- as.data.frame(pheno[,catVar_selected]);
		colnames(indSexCond) <- name_catVar_selected;
		indSexCond$labels <- rownames(pheno);
		rownames(indSexCond) <- rownames(pheno);

	};
	cond1 <- "Genotype" %in% colnames(indSexCond);
	cond2 <- "Sex:: Sex" %in% colnames(indSexCond);
	if (cond1==TRUE & cond2==TRUE) {
		indSexCond$cond_Geno_Sex <- paste0(indSexCond[,"Genotype"],"_",indSexCond[,"Sex:: Sex"]);
	};

	assign("indSexCond",indSexCond,.GlobalEnv);
};


#' Preformatting- random colour assignation
#'
#' @name random_colour_assignation.function 
#'
#' @param phenoColourDf Dataframe containing at least 2 columns, 
#' 1 "labels", containing the possible categories. Ex Genotype het, wt, homo
#' 2 "colPalette" hex colour associated to each varible. 
#' Ex for label="heterozygous", use colPalette="#8B0000" 
#' 
#' @return Several objects, the input for the full model phenoSocial_4analysis,
#' the imput for the selected model without highly corr variables, (noted as Sel1), 
#' phenoAll_4analysis_woHighCor_sel1. Both with and without Genotype.
#' The correlationMatrix, and other parameters for downstream functions.
#'
#' @examples
#' random_colour_assignation.function(phenoColourDf),
#'
#' @export
random_colour_assignation.function <- function(phenoColourDf){
	# description Function to assign a colour per variable for the plotting.

	#library("randomcoloR"); #give you the random colours in hex format
	nb <- length(unique(phenoColourDf[,1]));
	palette <- distinctColorPalette(nb);
	paletteRGBinput<- col2rgb(palette, alpha = FALSE);

	#changing the colours from hex format to RGB
	#palleteRGB <- NULL
	#for (case in 1:dim(paletteRGBinput)[2]){
	#	palleteRGB[case]<- paste0(paletteRGBinput[,case][1],",", paletteRGBinput[,case][2],",",paletteRGBinput[,case][3]);
	#	case=case+1;
	#};

	#phenoColourDf$colPalette <- palleteRGB;
	phenoColourDf$colPalette <- palette;
	assign("phenoColourDf",intramodelCombinationDfF ,.GlobalEnv);
};

