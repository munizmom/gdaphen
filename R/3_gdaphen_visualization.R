

#' plot0 - Cummulative variance plot over 10 components PCA
#'
#' @name cumVariance_Plot.function 
#' 
#' @param ResPCAmixData Input file is the output from MFA.function.
#'
#' @param nameOutput Short naming for the output of the plot.
#' 
#' @return The plot 0 is generated in the cumulativeVariance folder 
#'
#' @examples
#' cumVariance_Plot.function(res_PCAmixData_phenoAll_4analysis,"res_PCAmixData_phenoAll_4analysis")
#' cumVariance_Plot.function(res_PCAmixData_phenoAll_4analysis_woHighCor_sel1,"res_PCAmixData_phenoAll_4analysis_woHighCor_sel1")
#' cumVariance_Plot.function(res_PCAmixData_parameters_more_30perc_informative_Genotype,"res_PCAmixData_parameters_more_30perc_informative_Genotype")
#' cumVariance_Plot.function(res_PCAmixData_parameters_more_30perc_informative_woGenotype,"res_PCAmixData_parameters_more_30perc_informative_woGenotype")
#'
#' @export

####################################################
##   Plot A or 0: Cumulative variance and proportion  
##   variance explained by each dimension of the PCA
####################################################
cumVariance_Plot.function <- function(ResPCAmixData,nameOutput){
	f_cumVar <- "cumulativeVariance/"
	dir.create(file.path(wd,"/",f_pheno_discrimination,f_results,f_cumVar), showWarnings = F);


	datacumVar <- as.data.frame(ResPCAmixData$eig);
	datacumVar$dimensions <- gsub("dim ","",rownames(datacumVar));
	datacumVar$dimensions <- factor(datacumVar$dimensions, levels= seq(1:length(datacumVar$dimensions)));

	#selecting the first top 10 dimnsions (if they exist)
	dimsMax<- nrow(datacumVar);
	if (dimsMax>10){
		dimsTop <- 10;
	} else{
		dimsTop <- dimsMax;
	};

	datacumVar.sel10 <- datacumVar[c(1:dimsTop),]
	ylim.tmp <- max(datacumVar.sel10[,'Cumulative'])
	if (ylim.tmp<100){
		y.lim <-round(ylim.tmp,digits=0);
	} else{
		y.lim <-100;
	}
    #creating & toring the plot in a variable "cumVariancePlot.sel10"
	cumVariancePlot.sel10 <- ggplot(datacumVar.sel10, aes(x=dimensions)) +geom_col(aes(y=Proportion), fill="#3ac754")+ 
	scale_y_continuous("% variance",breaks=c(0,10,25,50,75,y.lim),limits=c(0,y.lim)) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size=7),axis.title=element_text(size=10,face="bold"),legend.position="none",plot.title = element_text(size=12))+
	scale_x_discrete("",expand = waiver()) +labs( title=paste("Variance explained by top10 dimensions" ), caption = "Source: Y.Herault team @IGBMC") +   
	guides(size=FALSE);
	cumVariancePlot.sel10 <- cumVariancePlot.sel10 + geom_line(aes(y = Cumulative), size = 1.5, color="#277335", group = 1)+ 
	geom_text(vjust = 1.3,hjust=0.7, size=2.5, aes(x=dimensions, y=Proportion, label=round(Proportion,digits=1),angle = 0)) 
	
    #saving the plot directly into a pdf into the defined folder
    pdf(paste0( wd,f_pheno_discrimination,f_results,f_cumVar,nameOutput,"_cumulativeVariance_perDimension.pdf"),height=3,width=5);
	print(cumVariancePlot.sel10);
	dev.off();
    
    #saving the numerical results in a excel file 
	write.xlsx(datacumVar, file = paste0(wd,f_pheno_discrimination,f_results,f_cumVar,nameOutput, "_cumulativeVariance_perDimension.xlsx"));

    #assigning the created variables in the global env for further use
	assign(paste0("cumVar_",nameOutput),datacumVar,.GlobalEnv);
	assign("datacumVar.sel10",datacumVar.sel10,.GlobalEnv);

};



#' plot1 - All variables correlation in PCA partial axes
#'
#' @name plot1_all_variables_correlation.function 
#'
#' @param name name for the plot or specific to the session
#'
#' @param catVarName Column name of the qualitative
#'  variable or categorical variable that is your dependant 
#'  variable, the one for which one you want to identify the most discriminative variables  
#'
#' @param nameExperiment The name of the experiment together with the number of variables we included. 
#' if is all model, allthe variables. Instead if is the selected model is less variables. 
#' ex: 
#' 1. nameExperiment <- paste0(name,"_", length(colnames(phenoSocial_4analysis)), " variables"); #all variables
#' 2. nameExperiment <- paste0(name,"_", length(colnames(phenoAll_4analysis_woHighCor_sel1)), " variables"); · without highly correlated
#' 3. nameExperiment <- paste0(name,"_", " more informative variables ",length(colnames(parameters_more_30perc_informative_Genotype)) ) #selected variables
#' 
#' @param margin maargin to add to the x and y axis to fit properly the plots limits 
#' 
#' @param ncolLegend Number of columns the legend should be split to, to format properly for publications.
#'
#' @param nameModel Short naming for the model you are using.
#' Ex: "fullMdl" for the model including all variables;"Sel1" for the model not including the
#' highly correlated variables, "Sel30" for the model including only the variables contributing
#' more thana 30%
#' 
#' @return The plot 1 is generated in the plot1 folder in three different kind of organizations
#' to select the better for publications (landscape, vertical, 1 plot alone or all plots together.
#'
#' @examples
#' plot1_all_variables_correlation.function(name,catVarName,nameExperiment,0.2,4,nameModel);
#'
#' @export
plot1_all_variables_correlation.function <- function(name,catVarName,nameExperiment,margin,ncolLegend,nameModel){
    load( file=paste0(wd, f_pheno_discrimination,f_results, f_Rdata,name,"_",nameModel,"_",catVarName,"_Robject_for_plots.RData")) #namesVariablesDf, df1, df2, df3,df4, df5, df6, MFA1, MFA2, MFA3,nbVarF, nbVarTotalF,
	df1$x <-0
	#df4$slope <-diff(c(0,df4[,2])) / diff(c(0,df4[,1])) #not needed
	Df1plot <- df1
	Df1plotlegend <- df1
	Df1plot$labels <- factor(Df1plot$labels, levels=unique(Df1plot$labels))
	Df1plot$colPalette <- factor(Df1plot$colPalette, levels=unique(Df1plot$colPalette))
	Df1plot$labels <- factor(Df1plot$labels, levels=unique(Df1plot$labels));
	Df1plot$colPalette <- factor(Df1plot$colPalette, levels=unique(Df1plot$colPalette));
	Df1plot$labels_nb <- factor(Df1plot$labels_nb, levels=unique(Df1plot$labels_nb))

	xylimp1=max(c(Df1plot[,1],Df1plot[,2]))+0.4+margin;
	mxylimp1=min(c(Df1plot[,1],Df1plot[,2]))-0.4-margin;

	p1 <- ggplot(Df1plot, aes(x=`dim 1`, y=`dim 2`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(mxylimp1,xylimp1))  + scale_y_continuous(limits=c(mxylimp1,xylimp1))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 2 (",MFA2, "% )"))+ theme(legend.position="none");
	p1 <- p1+ geom_segment(aes(x=x , y=x, xend=`dim 1`, yend=`dim 2`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=levels(Df1plot$colPalette),labels=levels(Df1plot$labels)) + ggforce::geom_circle(inherit.aes = FALSE, aes(x0=0, y0=0,r=xylimp1-0.4))+coord_fixed();
	     
	#plot1 label
	xylimp1Label=max(c(Df1plot[,1],Df1plot[,2]))+0.7+margin
	mxylimp1Label=min(c(Df1plot[,1],Df1plot[,2]))-0.7-margin

	p1label <- ggplot(Df1plot, aes(x=`dim 1`, y=`dim 2`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(mxylimp1Label,xylimp1Label))  + scale_y_continuous(limits=c(mxylimp1Label,xylimp1Label))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 2 (",MFA2, "% )"))+ theme(legend.position="none");
	p1label <- p1label+ geom_segment(aes(x=x , y=x, xend=`dim 1`, yend=`dim 2`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=levels(Df1plot$colPalette),labels=levels(Df1plot$labels)) + ggforce::geom_circle(inherit.aes = FALSE, aes(x0=0, y0=0,r=xylimp1Label-0.75))+coord_fixed();
	p1labelFull <- p1label + geom_text(aes(x=`dim 1`, y=`dim 2`,label=labels,size=0.1,vjust = 0.7,hjust=-0.13, angle=45))
	p1label_nb <- p1label + geom_text(aes(x=`dim 1`, y=`dim 2`,label=nb,size=0.1,vjust = ifelse(`dim 2`>0, -.2,`dim 2`+2.4),hjust=ifelse(`dim 2`>0, `dim 1`-0.1,`dim 1`+0.8), angle=45)) #ifelse(`dim 1`>0, 0.3,`dim 1`-0.5)

	#plot2
	xylimp2=max(c(Df1plot[,1],Df1plot[,3]))+0.3+margin;
	mxylimp2=min(c(Df1plot[,1],Df1plot[,3]))-0.3-margin;
	p2 <- ggplot(Df1plot, aes(x=`dim 1`, y=`dim 3`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(mxylimp2,xylimp2))  + scale_y_continuous(limits=c(mxylimp2,xylimp2))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none");
	p2 <- p2+ geom_segment(aes(x=x , y=x, xend=`dim 1`, yend=`dim 3`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=levels(Df1plot$colPalette),labels=levels(Df1plot$labels)) + ggforce::geom_circle(inherit.aes = FALSE, aes(x0=0, y0=0,r=xylimp2-0.4))+coord_fixed();
	#plot2 labzl
	xylimp2Label=max(c(Df1plot[,1],Df1plot[,3]))+0.7+margin;
	mxylimp2Label=min(c(Df1plot[,1],Df1plot[,3]))-0.7-margin;

	p2label <- ggplot(Df1plot, aes(x=`dim 1`, y=`dim 3`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(mxylimp2Label,xylimp2Label))  + scale_y_continuous(limits=c(mxylimp2Label,xylimp2Label))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none");
	p2label <- p2label+ geom_segment(aes(x=x , y=x, xend=`dim 1`, yend=`dim 3`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=levels(Df1plot$colPalette),labels=levels(Df1plot$labels)) + ggforce::geom_circle(inherit.aes = FALSE, aes(x0=0, y0=0,r=xylimp2Label-0.5))+coord_fixed();
	p2labelFull <- p2label + geom_text(aes(x=`dim 1`, y=`dim 3`,label=labels,size=0.1,vjust = 0.7,hjust=-0.13, angle=45));
	p2label_nb <- p2label + geom_text(aes(x=`dim 1`, y=`dim 3`,label=nb,size=0.1,vjust = ifelse(`dim 3`>0, -.2,`dim 3`+2.4),hjust=ifelse(`dim 3`>0, `dim 1`-0.1,`dim 1`+0.8), angle=45)) #ifelse(`dim 1`>0, 0.3,`dim 1`-0.5)

	##plot3
	xylimp3=max(c(Df1plot[,2],Df1plot[,3]))+0.5+margin;
	mxylimp3=min(c(Df1plot[,2],Df1plot[,3]))-0.5-margin;
	p3 <- ggplot(Df1plot, aes(x=`dim 2`, y=`dim 3`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(mxylimp3,xylimp3))  + scale_y_continuous(limits=c(mxylimp3,xylimp3))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 2 (",MFA2, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none");
	p3 <- p3+ geom_segment(aes(x=x , y=x, xend=`dim 2`, yend=`dim 3`), arrow=arrow(length=unit(0.3,"cm")))+ 
	scale_color_manual(values=levels(Df1plot$colPalette),labels=levels(Df1plot$labels))+ ggforce::geom_circle(inherit.aes = FALSE, aes(x0=0, y0=0,r=xylimp3-0.5))+coord_fixed();
	#plot" label
	xylimp3Label=max(c(Df1plot[,2],Df1plot[,3]))+0.8+margin;
	mxylimp3Label=min(c(Df1plot[,2],Df1plot[,3]))-0.8-margin;

	p3label <- ggplot(Df1plot, aes(x=`dim 2`, y=`dim 3`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(mxylimp3Label,xylimp3Label))  + scale_y_continuous(limits=c(mxylimp3Label,xylimp3Label))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 2 (",MFA2, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none");
	p3label <- p3label+ geom_segment(aes(x=x , y=x, xend=`dim 2`, yend=`dim 3`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=levels(Df1plot$colPalette),labels=levels(Df1plot$labels)) +ggforce::geom_circle(inherit.aes = FALSE, aes(x0=0, y0=0,r=xylimp3Label-0.7))+coord_fixed();
	p3labelFull <- p3label + geom_text(aes(x=`dim 2`, y=`dim 3`,label=labels,size=0.1,vjust = 0.7,hjust=-0.13, angle=45));
	p3label_nb <- p3label + geom_text(aes(x=`dim 2`, y=`dim 3`,label=nb,size=0.1,vjust = ifelse(`dim 3`>0, -.2,`dim 3`+2.4),hjust=ifelse(`dim 3`>0, `dim 2`-0.1,`dim 2`+0.8), angle=45)) #ifelse(`dim 1`>0, 0.3,`dim 1`-0.5)

	#legend
	pc_qp6<-qplot(data=Df1plotlegend, x=`dim 2`, y=`dim 3`,color=labels )+
		 	 geom_point(size=4,stat = "identity") + geom_hline(yintercept = 0, colour = "gray65") +
		 	  geom_vline(xintercept = 0, colour = "gray65")+ 
		 	  scale_color_manual(values=Df1plotlegend$colPalette,labels=Df1plotlegend$labels)+
		 	    labs(x = paste0("dim 2 (",MFA2, "% )"), 
		 	    	y=paste0("dim 3 (",MFA3, "% )")) +
		 	    guides(colour = guide_legend("labels",ncol = ncolLegend, byrow = FALSE))+ 
		 	    theme_gray(base_size = 12)+
		 	    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
		 	legend.title = element_text(size=12, face="bold"),legend.position="bottom", 
		 	legend.box = "horizontal",legend.text=element_text(size=10,face="bold"),
		 	legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line"));
	pc_qp6Vert<-qplot(data=Df1plotlegend, x=`dim 2`, y=`dim 3`,color=labels )+
		 	 geom_point(size=4,stat = "identity") + geom_hline(yintercept = 0, colour = "gray65") +
		 	  geom_vline(xintercept = 0, colour = "gray65")+ 
		 	  scale_color_manual(values=Df1plotlegend$colPalette,labels=Df1plotlegend$labels)+
		 	    labs(x = paste0("dim 2 (",MFA2, "% )"), 
		 	    	y=paste0("dim 3 (",MFA3, "% )")) +
		 	    guides(colour = guide_legend("labels",ncol = ncolLegend, byrow = FALSE))+ 
		 	    theme_gray(base_size = 12)+
		 	    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
		 	legend.title = element_text(size=12, face="bold"),legend.position="bottom", 
		 	legend.box = "vertical",legend.text=element_text(size=10,face="bold"),
		 	legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line"));


	pc_qp6_nb<-qplot(data=Df1plotlegend, x=`dim 2`, y=`dim 3`,color=labels_nb )+
		 	 geom_point(size=4,stat = "identity") + geom_hline(yintercept = 0, colour = "gray65") +
		 	  geom_vline(xintercept = 0, colour = "gray65")+ 
		 	  scale_color_manual(values=Df1plotlegend$colPalette,labels=Df1plotlegend$labels_nb)+
		 	    labs(x = paste0("dim 2 (",MFA2, "% )"), 
		 	    	y=paste0("dim 3 (",MFA3, "% )")) +
		 	    guides(colour = guide_legend("labels",ncol = ncolLegend, byrow = FALSE))+ 
		 	    theme_gray(base_size = 12)+
		 	    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
		 	legend.title = element_text(size=12, face="bold"),legend.position="bottom", 
		 	legend.box = "horizontal",legend.text=element_text(size=10,face="bold"),
		 	legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line"));
	#labels
	upText <-paste0("    ",nameExperiment,": input file with ",nbVarTotalF," phenotypic variables");
	upText.p <- ggparagraph(text = upText, face = "bold", size = 14, color = "black");
	midText <-paste0("    ",nameExperiment,": MFA considering ",nbVarF," non-correlated variables");
	midText.p <- ggparagraph(text = midText, face = "bold", size = 14, color = "black");

	pc_qp6_nb_vertical<-qplot(data=Df1plotlegend, x=`dim 2`, y=`dim 3`,color=labels_nb )+
		 	 geom_point(size=4,stat = "identity") + geom_hline(yintercept = 0, colour = "gray65") +
		 	  geom_vline(xintercept = 0, colour = "gray65")+ 
		 	  scale_color_manual(values=Df1plotlegend$colPalette,labels=Df1plotlegend$labels_nb)+
		 	    labs(x = paste0("dim 2 (",MFA2, "% )"), 
		 	    	y=paste0("dim 3 (",MFA3, "% )")) +
		 	    guides(colour = guide_legend("labels",ncol = 1, byrow = FALSE))+ 
		 	    theme_gray(base_size = 12)+
		 	    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
		 	legend.title = element_text(size=12, face="bold"),legend.position="bottom", 
		 	legend.box = "horizontal",legend.text=element_text(size=10,face="bold"),
		 	legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line"));
	#labels
	upText <-paste0("    ",nameExperiment,": input file with ",nbVarTotalF," phenotypic variables");
	upText.p <- ggparagraph(text = upText, face = "bold", size = 14, color = "black");
	midText <-paste0("    ",nameExperiment,": MFA considering ",nbVarF," non-correlated variables");
	midText.p <- ggparagraph(text = midText, face = "bold", size = 14, color = "black");

	# Arrange the plots on the same page  --> rats and all other plots x11(height=5, width=10.2)
	#X11(height=11.3, width=8.2) #8.3 x 11.7 i
	legend <- cowplot::get_legend(pc_qp6); 
	legendV <- cowplot::get_legend(pc_qp6Vert); 
	legend_nb <- cowplot::get_legend(pc_qp6_nb); 
	legend_vertical <- cowplot::get_legend(pc_qp6_nb_vertical); 
	blank<-rectGrob(gp=gpar(col="white")) # make a white spacer grob
	p1a <- ggarrange(upText.p, ncol = 1, nrow = 1);
	p1b <- ggarrange(midText.p, ncol = 1, nrow = 1);
	p2a <- ggarrange(p1, p2, p3,  ncol = 3, nrow = 1, heights = c(1, 1, 1));
	p12 <- ggarrange(p1a,p1b, p2a, ncol = 1, nrow = 3,heights = c(0.1,0.1, 1));
	p4b <- ggarrange(p1labelFull, p2labelFull, p3labelFull +  theme(legend.position="none"),  ncol = 3, nrow = 1, heights = c(1, 1, 1));
	p4_nb <- ggarrange(p1label_nb, p2label_nb, p3label_nb +  theme(legend.position="none"),  ncol = 3, nrow = 1, heights = c(1, 1, 1));

    p2v2b <-ggarrange(p1labelFull, p2labelFull, p3labelFull +  theme(legend.position="none"),  ncol = 1, nrow = 3, widths = c(1, 1, 1));
    p2v2c <-ggarrange(p1, p2, p3 +  theme(legend.position="none"),  ncol = 1, nrow = 3, widths = c(1, 1, 1));
    p2v2_nb <-ggarrange(p1label_nb, p2label_nb, p3label_nb +  theme(legend.position="none"),  ncol = 1, nrow = 3, widths = c(1.2, 1.2, 1.2));
    p2v2_nb_hor <-ggarrange(p1label_nb, p2label_nb, p3label_nb +  theme(legend.position="none"),  ncol = 3, nrow = 1, widths = c(1,1,1));
	p12v2b <- ggarrange(p1a,p1b,p2v2b,legendV, ncol = 1, nrow = 4,heights = c(0.02,0.02,1,0.5));
	p12v2c <- ggarrange(p1a,p1b,p2v2c,legendV, ncol = 1, nrow = 4,heights = c(0.02,0.02,1,0.5));
	p12v2_nb <- ggarrange(p1a,p1b,p2v2_nb,legend_nb, ncol = 1, nrow = 4,heights = c(0.02,0.02,1,0.5));
	p12v2_nb2 <- ggarrange(p1a,p1b,p2v2_nb,legend_vertical, ncol = 1, nrow = 4,heights = c(0.03,0.03,1.4,0.9));
	p12v2_nb_hor <- ggarrange(blank,p1a,p1b,p2v2_nb_hor,legend_nb, ncol = 1, nrow = 5,heights = c(0.1,0.1,0.1,1.2,1));

	#individual labelled plots
	p1labelarranged <- ggarrange(blank,p1a,p1b,p1labelFull, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p2labelarranged <- ggarrange(blank,p1a,p1b,p2labelFull, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p3labelarranged <- ggarrange(blank,p1a,p1b,p3labelFull, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));

	p1labelarranged_nb <- ggarrange(blank,p1a,p1b,p1label_nb, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p2labelarranged_nb <- ggarrange(blank,p1a,p1b,p2label_nb, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p3labelarranged_nb <- ggarrange(blank,p1a,p1b,p3label_nb, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));

	p1labelarranged_nb <- ggarrange(blank,p1a,p1b,p1label_nb, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p2labelarranged_nb <- ggarrange(blank,p1a,p1b,p2label_nb, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p3labelarranged_nb <- ggarrange(blank,p1a,p1b,p3label_nb, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));

	#final grouped objects
	pf1b <- ggarrange(blank,p12,blank,p4b,legend, ncol = 1,nrow = 5,heights = c(0.1,1,0.1,1,0.5)); #with/wo mabels
	pf2 <- ggarrange(blank,p12,legend, ncol = 1,nrow = 3,heights = c(0.1,1,0.5)); #just wo labels
	p4v2b <-	ggarrange(blank,p12v2b, ncol = 1,nrow = 2,heights = c(0.05,1)); #just with labels
	p4v2c <-	ggarrange(blank,p12v2c, ncol = 1,nrow = 2,heights = c(0.05,1)); #just with labels
	pf1_nb <- ggarrange(blank,p12,blank,p4_nb,legend_nb, ncol = 1,nrow = 5,heights = c(0.1,1,0.1,1,0.5)); #with/wo mabels
	p4v2_nb <-	ggarrange(blank,p12v2_nb, ncol = 1,nrow = 2,heights = c(0.05,1)); #just with labels
	p4v2_nb2 <-	ggarrange(blank,p12v2_nb2, ncol = 1,nrow = 2,heights = c(0.05,1)); #just with labels
	
	#saving into pdf final plots arranged
    
	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot1,nameExperiment, "_partialAxes_groupCor_IndividualPlot.pdf"),height=7, width=11);
	print(p1labelarranged);
	print(p2labelarranged);
	print(p3labelarranged);
	print(p1labelarranged_nb);
	print(p2labelarranged_nb);
	print(p3labelarranged_nb);
	dev.off();

	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot1,nameExperiment, "_partialAxes_groupCor_labelGrouped_verticalvlegend",".pdf"),height=23, width=12);
	print(p4v2_nb2)
	dev.off();

	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot1,nameExperiment, "_partialAxes_groupCor_labelGrouped_vertical",".pdf"),height=20, width=12);
	print(p4v2_nb2)
	print(p4v2_nb);
	print(p4v2b);
	print(p4v2c);
	dev.off();
	
	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot1,nameExperiment, "_partialAxes_groupCor_labelGrouped_horizontal",".pdf"),height=8, width=14);
	print(p12v2_nb_hor)
	print(pf2);
	dev.off();

	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot1,nameExperiment, "_partialAxes_groupCor_allGrouped",".pdf"),height=18, width=20);
	print(pf1_nb);
	print(pf1b);
	dev.off();


	};




#' plot2 - Grouped variables correlation in PCA typical plot
#'
#' @name plot2_grouped_vars_correlation_PCAs.function 
#'
#' @param name name for the plot or specific to the session
#'
#' @param catVarName Column name of the qualitative
#'  variable or categorical variable that is your dependant 
#'  variable, the one for which one you want to identify the most discriminative variables  
#'
#' @param nameExperiment The name of the experiment together with the number of variables we included. 
#' if is all model, allthe variables. Instead if is the selected model is less variables. 
#' ex: 
#' 1. nameExperiment <- paste0(name,"_", length(colnames(phenoSocial_4analysis)), " variables"); #all variables
#' 2. nameExperiment <- paste0(name,"_", length(colnames(phenoAll_4analysis_woHighCor_sel1)), " variables"); · without highly correlated
#' 3. nameExperiment <- paste0(name,"_", " more informative variables ",length(colnames(parameters_more_30perc_informative_Genotype)) ) #selected variables
#' 
#' @param nrowLegend Number of rows the legend should be split to, to format properly for publications.
#'
#' @param nameModel Short naming for the model you are using.
#' Ex: "fullMdl" for the model including all variables;"Sel1" for the model not including the
#' highly correlated variables, "Sel30" for the model including only the variables contributing
#' more thana 30%
#' 
#' @return The plot 2 is generated in the plot2 folder in three different kind of organizations
#' to select the better for publications (landscape, vertical, 1 plot alone or all plots together.
#'
#' @examples
#' plot2_grouped_vars_correlation_PCAs.function(name,catVarName,nameExperiment,4,nameModel);
#'
#' @export
plot2_grouped_vars_correlation_PCAs.function <- function(name,catVarName,nameExperiment,nrowLegend,nameModel){
    load( file=paste0(wd, f_pheno_discrimination,f_results, f_Rdata,name,"_",nameModel,"_",catVarName,"_Robject_for_plots.RData")) #namesVariablesDf, df1, df2, df3,df4, df5, df6, MFA1, MFA2, MFA3,nbVarF, nbVarTotalF,
	Df2plot <- df2
	Df2plot$labels <- factor(Df2plot$labels, levels=unique(Df2plot$labels))

	#x and y axis limits
	min_value_3d <- min(c(Df2plot[,1],Df2plot[,2],Df2plot[,3])) 
	max_value_3d <- max(c(Df2plot[,1],Df2plot[,2],Df2plot[,3])) 

	#2) plots for the variables selected on the analysis
	pc_qp1<- qplot(x=`dim 1`, y=`dim 2`, data=as.data.frame(Df2plot), color=labels)+  scale_x_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.02))  + scale_y_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.02)) + theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ scale_color_manual(breaks = labels, values=Df2plot$colPalette)+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 2 (",MFA2, "% )"))+ theme(legend.position="none");  #take out the ,labels=my_leg1 to not put the italics
	pc_qp1label<-pc_qp1 +geom_text(aes(label=labels,size=2,vjust = -0.75,hjust=0.4,angle=10))   #vjust = 0.1,

	pc_qp2<-qplot(x=`dim 1`, y=`dim 3`, data=as.data.frame(Df2plot), color=labels) + scale_x_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.02))  + scale_y_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.02)) + theme_gray(base_size = 11)+ theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ scale_color_manual(breaks = labels, values=Df2plot$colPalette)+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none"); 
	pc_qp2label<-pc_qp2 +geom_text(aes(label=labels,size=2,vjust = -0.75,hjust=0.4,angle=10)) 

	pc_qp3<-qplot(x=`dim 2`, y=`dim 3`, data=as.data.frame(Df2plot), color=labels) + scale_x_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.02))  + scale_y_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.02)) + theme_gray(base_size = 11)+ theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ scale_color_manual(breaks = labels, values=Df2plot$colPalette)+  labs(x = paste0("dim 2 (",MFA2, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none"); 
	pc_qp3label<-pc_qp3 +geom_text(aes(label=labels,size=2,vjust = -0.75,hjust=0.4,angle=10)) 
	pc_qp6<-qplot(data=Df2plot, x=`dim 2`, y=`dim 3`,color=labels )+
		 	 geom_point(size=4,stat = "identity") + geom_hline(yintercept = 0, colour = "gray65") +
		 	  geom_vline(xintercept = 0, colour = "gray65")+ 
		 	  scale_color_manual(values=Df2plot$colPalette,labels=Df2plot$labels)+
		 	    labs(x = paste0("dim 2 (",MFA2, "% )"), 
		 	    	y=paste0("dim 3 (",MFA3, "% )")) +
		 	    guides(colour = guide_legend("labels",nrow = nrowLegend, byrow = TRUE))+ 
		 	    theme_gray(base_size = 12)+
		 	    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
		 	legend.title = element_text(size=12, face="bold"),legend.position="bottom", 
		 	legend.box = "horizontal",legend.text=element_text(size=10,face="bold"),
		 	legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line"));

	#labels
	upText <-paste0("    ",nameExperiment,": input file with ",nbVarTotalF," phenotypic variables");
	upText.p <- ggparagraph(text = upText, face = "bold", size = 14, color = "black");
	midText <-paste0("    ",nameExperiment,": MFA considering ",nbVarF," non-correlated variables");
	midText.p <- ggparagraph(text = midText, face = "bold", size = 14, color = "black");


	# Arrange the plots on the same page  --> rats and all other plots x11(height=5, width=10.2)
	#X11(height=13.3, width=20.2) #8.3 x 11.7 
	legend <- cowplot::get_legend(pc_qp6); 
	blank<-rectGrob(gp=gpar(col="white")) # make a white spacer grob
	p1a <- ggarrange(upText.p, ncol = 1, nrow = 1);
	p1b <- ggarrange(midText.p, ncol = 1, nrow = 1);
	p2 <- ggarrange(pc_qp1, pc_qp2, pc_qp3,  ncol = 3, nrow = 1, heights = c(1, 1, 1));
	p12 <- ggarrange(p1a,p1b, p2, ncol = 1, nrow = 3,heights = c(0.1,0.1, 1));
	p4 <- ggarrange(pc_qp1label, pc_qp2label, pc_qp3label +  theme(legend.position="none"),  ncol = 3, nrow = 1, heights = c(1, 1, 1));


    p2v2 <-ggarrange(pc_qp1label, pc_qp2label, pc_qp3label +  theme(legend.position="none"),  ncol = 1, nrow = 3, widths = c(1, 1, 1));
	p12v2 <- ggarrange(p1a,p1b,p2v2, ncol = 1, nrow = 3,heights = c(0.02,0.02,1));
	#p4v2 <-	ggarrange(blank,p12v2,legend, ncol = 1,nrow = 3,heights = c(0.1,1,0.5));w legend
	#p4v2 <-	ggarrange(blank,p12v2, ncol = 1,nrow = 2,heights = c(0.05,1)); #wo legend

	#individual labelled plots
	pc_qp1labelarranged <- ggarrange(blank,p1a,p1b,pc_qp1label, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	pc_qp2labelarranged <- ggarrange(blank,p1a,p1b,pc_qp2label, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	pc_qp3labelarranged <- ggarrange(blank,p1a,p1b,pc_qp3label, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));

	#final grouped objects
	pf1 <- ggarrange(blank,p12,blank,p4,legend, ncol = 1,nrow = 5,heights = c(0.1,1,0.1,1,0.5)); #with/wo mabels
	pf2 <- ggarrange(blank,p12,legend, ncol = 1,nrow = 3,heights = c(0.1,1,0.5)); #just wo labels
	p4v2 <-	ggarrange(blank,p12v2, ncol = 1,nrow = 2,heights = c(0.05,1)); #just with labels

	#saving into pdf final plots arranged
	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot2,nameExperiment, "_GroupsContribution_IndividualPlot.pdf"),height=7, width=11)
	print(pc_qp1labelarranged);
	print(pc_qp2labelarranged);
	print(pc_qp3labelarranged);
	dev.off();

	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot2,nameExperiment, "_GroupsContribution_labelGrouped",".pdf"),height=15, width=12)
	print(p4v2);
	dev.off();
	
	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot2,nameExperiment, "_GroupsContribution_unlabelGrouped",".pdf"),height=6, width=14)
	print(pf2);
	dev.off();

	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot2,nameExperiment, "_GroupsContribution_allGrouped",".pdf"),height=14, width=30)
	print(pf1);
	dev.off();

};





#' plot3 - Individuals observations in PCA
#'
#' @name plot3_Individual_component_PCA.function 
#'
#' @param name name for the plot or specific to the session
#'
#' @param catVarName Column name of the qualitative
#'  variable or categorical variable that is your dependant 
#'  variable, the one for which one you want to identify the most discriminative variables  
#' 
#' @param nameExperiment The name of the experiment together with the number of variables we included. 
#' if is all model, allthe variables. Instead if is the selected model is less variables. 
#' ex: 
#' 1. nameExperiment <- paste0(name,"_", length(colnames(phenoSocial_4analysis)), " variables"); #all variables
#' 2. nameExperiment <- paste0(name,"_", length(colnames(phenoAll_4analysis_woHighCor_sel1)), " variables"); · without highly correlated
#' 3. nameExperiment <- paste0(name,"_", " more informative variables ",length(colnames(parameters_more_30perc_informative_Genotype)) ) #selected variables
#'
#' @param condition4Plot If several conditions can be plotted. Ex Genotype, Genotype+ treatments, etc.
#' identify here the name of the variable to plot, is the name of the column from df3 that you want to plot.
#'
#' @param colPalette Column name from df3 that have the color assignation to plot the columns condition4Plot
#'
#' @param nrowLegend Number of rows the legend should be split to, to format properly for publications.
#'
#' @param nameModel Short naming for the model you are using.
#' Ex: "fullMdl" for the model including all variables;"Sel1" for the model not including the
#' highly correlated variables, "Sel30" for the model including only the variables contributing
#' more thana 30%
#' 
#' @return The plot 3 is generated in the plot3 folder in three different kind of organizations
#' to select the better for publications (landscape, vertical, 1 plot alone or all plots together.
#'
#' @examples
#' plot3_Individual_component_PCA.function(name,catVarName,"Genotype","colPalette_cond_Geno",nameExperiment,4,nameModel) #run for all the genotypes comparison
#' plot3_Individual_component_PCA.function(name,catVarName,"cond_Geno_Sex","colPalette_cond_Geno_Sex",paste0(nameExperiment,"_condGeno_sex"),4,nameModel) #run for all the genotypes comparison
#'
#' @export
plot3_Individual_component_PCA.function <- function(name,catVarName,condition4Plot,colPalette,nameExperiment,nrowLegend,nameModel){
    load( file=paste0(wd, f_pheno_discrimination,f_results, f_Rdata,name,"_",nameModel,"_",catVarName,"_Robject_for_plots.RData")) #namesVariablesDf, df1, df2, df3,df4, df5, df6, MFA1, MFA2, MFA3,nbVarF, nbVarTotalF,

    Df3plot <- df3
	Df3plot$labels <- factor(Df3plot$labels, levels=unique(Df3plot$labels))
	Df3plot$conditionToPlot <- Df3plot[,condition4Plot]
	Df3plot$colPalette <-Df3plot[,colPalette]
	#now creating new duplciated fake columns to be able to plot succesfully in ggplot using the default ggplot col names
	Df3plot$conditionToPlot <- factor(Df3plot$conditionToPlot, levels=unique(Df3plot$conditionToPlot))

	#x and y axis limits
	min_value_3d <- min(c(Df3plot[,1],Df3plot[,2],Df3plot[,3]))
	max_value_3d <- max(c(Df3plot[,1],Df3plot[,2],Df3plot[,3])) 

	#2) plots for the variables selected on the analysis
	pc_qp1<- qplot(x=`dim 1`, y=`dim 2`, data=as.data.frame(Df3plot), color=conditionToPlot) + scale_x_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.02))  + scale_y_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.02)) + theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ scale_color_manual(breaks = unique(Df3plot$conditionToPlot), values=unique(Df3plot$colPalette))+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 2 (",MFA2, "% )"))+ theme(legend.position="none");  #take out the ,labels=my_leg1 to not put the italics
	pc_qp1label<-pc_qp1 + geom_text(aes(label=labels,size=0.7,vjust = -0.75,hjust=0.4,angle=10)); #vjust = 0.1,


	pc_qp2<-qplot(x=`dim 1`, y=`dim 3`, data=as.data.frame(Df3plot), color=conditionToPlot) + scale_x_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.02))  + scale_y_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.02))  + theme_gray(base_size = 11)+ theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ scale_color_manual(breaks = unique(Df3plot$conditionToPlot), values=unique(Df3plot$colPalette))+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none"); 
	pc_qp2label <-pc_qp2 + geom_text(aes(label=labels,size=0.7,vjust = -0.75,hjust=0.4,angle=10));  #vjust = 0.1,


	pc_qp3<-qplot(x=`dim 2`, y=`dim 3`, data=as.data.frame(Df3plot), color=conditionToPlot) + scale_x_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.02))  + scale_y_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.02))  + theme_gray(base_size = 11)+ theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ scale_color_manual(breaks = unique(Df3plot$conditionToPlot), values=unique(Df3plot$colPalette))+  labs(x = paste0("dim 2 (",MFA2, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none"); 
	pc_qp3label<-pc_qp3 + geom_text(aes(label=labels,size=0.7,vjust = -0.75,hjust=0.4,angle=10));  #vjust = 0.1,

	pc_qp6<-qplot(data=Df3plot, x=`dim 2`, y=`dim 3`,color=conditionToPlot )+
		 	 geom_point(size=4,stat = "identity") + geom_hline(yintercept = 0, colour = "gray65") +
		 	  geom_vline(xintercept = 0, colour = "gray65")+ 
		 	  scale_color_manual(values=unique(Df3plot$colPalette),labels=levels(Df3plot$conditionToPlot))+
		 	    labs(x = paste0("dim 2 (",MFA2, "% )"), 
		 	    	y=paste0("dim 3 (",MFA3, "% )")) +
		 	    guides(colour = guide_legend("labels",nrow = nrowLegend, byrow = TRUE))+ 
		 	    theme_gray(base_size = 12)+
		 	    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
		 	legend.title = element_text(size=12, face="bold"),legend.position="bottom", 
		 	legend.box = "horizontal",legend.text=element_text(size=10,face="bold"),
		 	legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line"));

	#labelsf
	upText <-paste0("    ",nameExperiment,": input file with ",nbVarTotalF," phenotypic variables");
	upText.p <- ggparagraph(text = upText, face = "bold", size = 14, color = "black");
	midText <-paste0("    ",nameExperiment,": MFA considering ",nbVarF," non-correlated variables");
	midText.p <- ggparagraph(text = midText, face = "bold", size = 14, color = "black");

	# Arrange the plots on the same page  --> rats and all other plots x11(height=5, width=10.2)
	#X11(height=11.3, width=8.2) #8.3 x 11.7 i
	legend <- cowplot::get_legend(pc_qp6); 
	blank<-rectGrob(gp=gpar(col="white")) # make a white spacer grob
	p1a <- ggarrange(upText.p, ncol = 1, nrow = 1);
	p1b <- ggarrange(midText.p, ncol = 1, nrow = 1);
	p2 <- ggarrange(pc_qp1, pc_qp2, pc_qp3,  ncol = 3, nrow = 1, heights = c(1, 1, 1));
	p12 <- ggarrange(p1a,p1b, p2, ncol = 1, nrow = 3,heights = c(0.1,0.1, 1));
	p4 <- ggarrange(pc_qp1label, pc_qp2label, pc_qp3label +  theme(legend.position="none"),  ncol = 3, nrow = 1, heights = c(1, 1, 1));

    p2v2 <-ggarrange(pc_qp1label, pc_qp2label, pc_qp3label +  theme(legend.position="none"),  ncol = 1, nrow = 3, widths = c(1, 1, 1));
	p12v2 <- ggarrange(p1a,p1b,p2v2, ncol = 1, nrow = 3,heights = c(0.02,0.02,1));
	#p4v2 <-	ggarrange(blank,p12v2,legend, ncol = 1,nrow = 3,heights = c(0.1,1,0.5));w legend
	#p4v2 <-	ggarrange(blank,p12v2, ncol = 1,nrow = 2,heights = c(0.05,1)); #wo legend

	#individual labelled plots
	pc_qp1labelarranged <- ggarrange(blank,p1a,p1b,pc_qp1label, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	pc_qp2labelarranged <- ggarrange(blank,p1a,p1b,pc_qp2label, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	pc_qp3labelarranged <- ggarrange(blank,p1a,p1b,pc_qp3label, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));

	#final grouped objects
	pf1 <- ggarrange(blank,p12,blank,p4,legend, ncol = 1,nrow = 5,heights = c(0.1,1,0.1,1,0.5)); #with/wo mabels
	pf2 <- ggarrange(blank,p12,legend, ncol = 1,nrow = 3,heights = c(0.1,1,0.5)); #just wo labels
	p4v2 <-	ggarrange(blank,p12v2, ncol = 1,nrow = 2,heights = c(0.05,1)); #just with labels

	#saving into pdf final plots arranged
	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot3,nameExperiment, "_individualCoordinates_IndividualPlot.pdf"),height=7, width=11)
	print(pc_qp1labelarranged);
	print(pc_qp2labelarranged);
	print(pc_qp3labelarranged);
	dev.off();

	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot3,nameExperiment, "_individualCoordinates_labelGrouped",".pdf"),height=20, width=12)
	print(p4v2);
	dev.off();
	
	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot3,nameExperiment, "_individualCoordinates_unlabelGrouped",".pdf"),height=8, width=14)
	print(pf2);
	dev.off();

	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot3,nameExperiment, "_individualCoordinates_allGrouped",".pdf"),height=18, width=30)
	print(pf1);
	dev.off();
};


#' plot4 - Squared loadings of the variables in PCA
#'
#' @name plot4_squared_loads_all_vars_PCA.function 
#'
#' @param name name for the plot or specific to the session
#'
#' @param catVarName Column name of the qualitative
#'  variable or categorical variable that is your dependant 
#'  variable, the one for which one you want to identify the most discriminative variables  
#'
#' @param nameExperiment The name of the experiment together with the number of variables we included. 
#' if is all model, allthe variables. Instead if is the selected model is less variables. 
#' ex: 
#' 1. nameExperiment <- paste0(name,"_", length(colnames(phenoSocial_4analysis)), " variables"); #all variables
#' 2. nameExperiment <- paste0(name,"_", length(colnames(phenoAll_4analysis_woHighCor_sel1)), " variables"); · without highly correlated
#' 3. nameExperiment <- paste0(name,"_", " more informative variables ",length(colnames(parameters_more_30perc_informative_Genotype)) ) #selected variables
#' 
#' @param nrowLegend Number of rows the legend should be split to, to format properly for publications.
#'
#' @param nameModel Short naming for the model you are using.
#' Ex: "fullMdl" for the model including all variables;"Sel1" for the model not including the
#' highly correlated variables, "Sel30" for the model including only the variables contributing
#' more thana 30%
#' 
#' @return The plot 4 is generated in the plot4 folder in three different kind of organizations
#' to select the better for publications (landscape, vertical, 1 plot alone or all plots together.
#'
#' @examples
#' plot4_squared_loads_all_vars_PCA.function(name,catVarName,nameExperiment,4,nameModel);
#'
#' @export
plot4_squared_loads_all_vars_PCA.function <- function(name,catVarName,nameExperiment,nrowLegend,nameModel){
    load( file=paste0(wd, f_pheno_discrimination,f_results, f_Rdata,name,"_",nameModel,"_",catVarName,"_Robject_for_plots.RData")) #namesVariablesDf, df1, df2, df3,df4, df5, df6, MFA1, MFA2, MFA3,nbVarF, nbVarTotalF,
	df4$x <-0
	Df4plot <- df4
	#Df4plot$labels <- factor(Df4plot$labels, levels=unique(Df4plot$labels))
	#Df4plot$colPalette <- factor(Df4plot$colPalette, levels=unique(Df4plot$colPalette))

	p1 <- ggplot(Df4plot, aes(x=`dim 1`, y=`dim 2`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(min(Df4plot[,1])-0.02,max(Df4plot[,1])+0.1))  + scale_y_continuous(limits=c(min(Df4plot[,2])-0.02,max(Df4plot[,2])+0.05))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 2 (",MFA2, "% )"))+ theme(legend.position="none");
	p1 <- p1+ geom_segment(aes(x=x , y=x, xend=`dim 1`, yend=`dim 2`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=Df4plot$colPalette,labels=Df4plot$labels) 

	p1label <- p1 + geom_text(aes(x=`dim 1`, y=`dim 2`,label=labels,size=0.7,vjust = 0.6,hjust=-0.01, angle=45))


	p2 <- ggplot(Df4plot, aes(x=`dim 1`, y=`dim 3`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(min(Df4plot[,1])-0.02,max(Df4plot[,1])+0.1))  + scale_y_continuous(limits=c(min(Df4plot[,3])-0.02,max(Df4plot[,3])+0.05))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none");
	p2 <- p2+ geom_segment(aes(x=x , y=x, xend=`dim 1`, yend=`dim 3`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=Df4plot$colPalette,labels=Df4plot$labels) 

	p2label <- p2 + geom_text(aes(x=`dim 1`, y=`dim 3`,label=labels,size=0.7,vjust = 0.6,hjust=-0.01, angle=45))

	p3 <- ggplot(Df4plot, aes(x=`dim 2`, y=`dim 3`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(min(Df4plot[,2])-0.02,max(Df4plot[,2])+0.1))  + scale_y_continuous(limits=c(min(Df4plot[,3])-0.02,max(Df4plot[,3])+0.05))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 2 (",MFA2, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none");
	p3 <- p3+ geom_segment(aes(x=x , y=x, xend=`dim 2`, yend=`dim 3`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=Df4plot$colPalette,labels=Df4plot$labels) 

	p3label <- p3 + geom_text(aes(x=`dim 2`, y=`dim 3`,label=labels,size=0.7,vjust = 0.6,hjust=-0.01, angle=45))

	#legend
	pc_qp6<-qplot(data=Df4plot, x=`dim 2`, y=`dim 3`,color=labels )+
		 	 geom_point(size=4,stat = "identity") + geom_hline(yintercept = 0, colour = "gray65") +
		 	  geom_vline(xintercept = 0, colour = "gray65")+ 
		 	  scale_color_manual(values=Df4plot$colPalette,labels=Df4plot$labels)+
		 	    labs(x = paste0("dim 2 (",MFA2, "% )"), 
		 	    	y=paste0("dim 3 (",MFA3, "% )")) +
		 	    guides(colour = guide_legend("labels",nrow = nrowLegend, byrow = TRUE))+ 
		 	    theme_gray(base_size = 12)+
		 	    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
		 	legend.title = element_text(size=12, face="bold"),legend.position="bottom", 
		 	legend.box = "horizontal",legend.text=element_text(size=10,face="bold"),
		 	legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line"));

	#labels
	upText <-paste0("    ",nameExperiment,": input file with ",nbVarTotalF," phenotypic variables");
	upText.p <- ggparagraph(text = upText, face = "bold", size = 14, color = "black");
	midText <-paste0("    ",nameExperiment,": MFA considering ",nbVarF," non-correlated variables");
	midText.p <- ggparagraph(text = midText, face = "bold", size = 14, color = "black");

	# Arrange the plots on the same page  --> rats and all other plots x11(height=5, width=10.2)
	#X11(height=11.3, width=8.2) #8.3 x 11.7 i
	legend <- cowplot::get_legend(pc_qp6); 
	blank<-rectGrob(gp=gpar(col="white")) # make a white spacer grob
	p1a <- ggarrange(upText.p, ncol = 1, nrow = 1);
	p1b <- ggarrange(midText.p, ncol = 1, nrow = 1);
	p2 <- ggarrange(p1, p2, p3,  ncol = 3, nrow = 1, heights = c(1, 1, 1));
	p12 <- ggarrange(p1a,p1b, p2, ncol = 1, nrow = 3,heights = c(0.1,0.1, 1));
	p4 <- ggarrange(p1label, p2label, p3label +  theme(legend.position="none"),  ncol = 3, nrow = 1, heights = c(1, 1, 1));

    p2v2 <-ggarrange(p1label, p2label, p3label +  theme(legend.position="none"),  ncol = 1, nrow = 3, widths = c(1, 1, 1));
	p12v2 <- ggarrange(p1a,p1b,p2v2, ncol = 1, nrow = 3,heights = c(0.02,0.02,1));

	#individual labelled plots
	p1labelarranged <- ggarrange(blank,p1a,p1b,p1label, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p2labelarranged <- ggarrange(blank,p1a,p1b,p2label, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p3labelarranged <- ggarrange(blank,p1a,p1b,p3label, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));

	#final grouped objects
	pf1 <- ggarrange(blank,p12,blank,p4,legend, ncol = 1,nrow = 5,heights = c(0.1,1,0.1,1,0.5)); #with/wo mabels
	pf2 <- ggarrange(blank,p12,legend, ncol = 1,nrow = 3,heights = c(0.1,1,0.5)); #just wo labels
	p4v2 <-	ggarrange(blank,p12v2, ncol = 1,nrow = 2,heights = c(0.05,1)); #just with labels

	#saving into pdf final plots arranged
	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot4,nameExperiment, "_sqload_IndividualPlot.pdf"),height=7, width=11)
	print(p1labelarranged);
	print(p2labelarranged);
	print(p3labelarranged);
	dev.off();

	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot4,nameExperiment, "_sqload_labelGrouped",".pdf"),height=20, width=9)
	print(p4v2);
	dev.off();
	
	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot4,nameExperiment, "_sqload_unlabelGrouped",".pdf"),height=8, width=14)
	print(pf2);
	dev.off();

	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot4,nameExperiment, "_sqload_allGrouped",".pdf"),height=18, width=20)
	print(pf1);
	dev.off();

};


#' plot5 -	Grouped variables correlation in PCA partial axes circle
#'
#' @name plot5_circle_correlation_all_vars.function 
#'
#' @param name name for the plot or specific to the session
#'
#' @param catVarName Column name of the qualitative
#'  variable or categorical variable that is your dependant 
#'  variable, the one for which one you want to identify the most discriminative variables  
#' 
#' @param nameExperiment The name of the experiment together with the number of variables we included. 
#' if is all model, allthe variables. Instead if is the selected model is less variables. 
#' ex: 
#' 1. nameExperiment <- paste0(name,"_", length(colnames(phenoSocial_4analysis)), " variables"); #all variables
#' 2. nameExperiment <- paste0(name,"_", length(colnames(phenoAll_4analysis_woHighCor_sel1)), " variables"); · without highly correlated
#' 3. nameExperiment <- paste0(name,"_", " more informative variables ",length(colnames(parameters_more_30perc_informative_Genotype)) ) #selected variables
#'
#' @param margin maargin to add to the x and y axis to fit properly the plots limits 
#' 
#' @param nrowLegend Number of rows the legend should be split to, to format properly for publications.
#'
#' @param nameModel Short naming for the model you are using.
#' Ex: "fullMdl" for the model including all variables;"Sel1" for the model not including the
#' highly correlated variables, "Sel30" for the model including only the variables contributing
#' more thana 30%
#' 
#' @return The plot 5 is generated in the plot5 folder in three different kind of organizations
#' to select the better for publications (landscape, vertical, 1 plot alone or all plots together.
#'
#' @examples
#' plot5_circle_correlation_all_vars.function(name,catVarName,nameExperiment,0.2,4,nameModel);
#'
#' @export
plot5_circle_correlation_all_vars.function <- function(name,catVarName,nameExperiment,margin,nrowLegend,nameModel){
    load( file=paste0(wd, f_pheno_discrimination,f_results, f_Rdata,name,"_",nameModel,"_",catVarName,"_Robject_for_plots.RData")) #namesVariablesDf, df1, df2, df3,df4, df5, df6, MFA1, MFA2, MFA3,nbVarF, nbVarTotalF,

    df5$x <-0
	Df5plot <- df5

	xylimp1=max(c(Df5plot[,1],Df5plot[,2]))+0.3 + margin;
	mxylimp1=min(c(Df5plot[,1],Df5plot[,2]))-0.3 - margin;
	p1 <- ggplot(Df5plot, aes(x=`dim 1`, y=`dim 2`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(mxylimp1,xylimp1))  + scale_y_continuous(limits=c(mxylimp1,xylimp1))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 2 (",MFA2, "% )"))+ theme(legend.position="none");
	p1 <- p1+ geom_segment(aes(x=x , y=x, xend=`dim 1`, yend=`dim 2`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=Df5plot$colPalette,labels=Df5plot$labels) +ggforce::geom_circle(inherit.aes = FALSE, aes(x0=0, y0=0,r=xylimp1-0.4))+coord_fixed()
	#plot1 labl
	xylimp1Label=max(c(Df5plot[,1],Df5plot[,2]))+0.7 + margin
	mxylimp1Label=min(c(Df5plot[,1],Df5plot[,2]))-0.7 - margin

	p1label <- ggplot(Df5plot, aes(x=`dim 1`, y=`dim 2`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(mxylimp1Label,xylimp1Label))  + scale_y_continuous(limits=c(mxylimp1Label,xylimp1Label))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 2 (",MFA2, "% )"))+ theme(legend.position="none");
	p1label <- p1label+ geom_segment(aes(x=x , y=x, xend=`dim 1`, yend=`dim 2`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=Df5plot$colPalette,labels=Df5plot$labels) +ggforce::geom_circle(inherit.aes = FALSE, aes(x0=0, y0=0,r=xylimp1Label-0.15))+coord_fixed()
	p1labelShort <- p1label + geom_text(aes(x=`dim 1`, y=`dim 2`,label=labelsParameter,size=0.1,vjust = 0.7,hjust=-0.13, angle=45))
	p1labelFull <- p1label + geom_text(aes(x=`dim 1`, y=`dim 2`,label=labels,size=0.1,vjust = 0.7,hjust=-0.13, angle=45))
	p1label_nb <- p1label + geom_text(aes(x=`dim 1`, y=`dim 2`,label=nb,size=0.1,vjust = ifelse(`dim 2`>0, -.2,`dim 2`+2.4),hjust=ifelse(`dim 2`>0, `dim 1`-0.1,`dim 1`+0.8), angle=45)) #ifelse(`dim 1`>0, 0.3,`dim 1`-0.5)

	#plot2
	xylimp2=max(c(Df5plot[,1],Df5plot[,3]))+0.3 + margin
	mxylimp2=min(c(Df5plot[,1],Df5plot[,3]))-0.3 - margin
	p2 <- ggplot(Df5plot, aes(x=`dim 1`, y=`dim 3`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(mxylimp2,xylimp2))  + scale_y_continuous(limits=c(mxylimp2,xylimp2))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none");
	p2 <- p2+ geom_segment(aes(x=x , y=x, xend=`dim 1`, yend=`dim 3`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=Df5plot$colPalette,labels=Df5plot$labels) +ggforce::geom_circle(inherit.aes = FALSE, aes(x0=0, y0=0,r=xylimp2-0.4))+coord_fixed()
	#plot2 labzl
	xylimp2Label=max(c(Df5plot[,1],Df5plot[,3]))+0.7 + margin
	mxylimp2Label=min(c(Df5plot[,1],Df5plot[,3]))-0.7 - margin

	p2label <- ggplot(Df5plot, aes(x=`dim 1`, y=`dim 3`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(mxylimp2Label,xylimp2Label))  + scale_y_continuous(limits=c(mxylimp2Label,xylimp2Label))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none");
	p2label <- p2label+ geom_segment(aes(x=x , y=x, xend=`dim 1`, yend=`dim 3`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=Df5plot$colPalette,labels=Df5plot$labels) +ggforce::geom_circle(inherit.aes = FALSE, aes(x0=0, y0=0,r=xylimp2Label-0.15))+coord_fixed()
	p2labelShort <- p2label + geom_text(aes(x=`dim 1`, y=`dim 3`,label=labelsParameter,size=0.1,vjust = 0.7,hjust=-0.13, angle=45))
	p2labelFull <- p2label + geom_text(aes(x=`dim 1`, y=`dim 3`,label=labels,size=0.1,vjust = 0.7,hjust=-0.13, angle=45))
	p2label_nb <- p2label + geom_text(aes(x=`dim 1`, y=`dim 3`,label=nb,size=0.1,vjust = ifelse(`dim 3`>0, -.2,`dim 3`+2.4),hjust=ifelse(`dim 3`>0, `dim 1`-0.1,`dim 1`+0.8), angle=45)) #ifelse(`dim 1`>0, 0.3,`dim 1`-0.5)
	
	##plot3
	xylimp3=max(c(Df5plot[,2],Df5plot[,3]))+0.3 + margin
	mxylimp3=min(c(Df5plot[,2],Df5plot[,3]))-0.3 - margin
	p3 <- ggplot(Df5plot, aes(x=`dim 2`, y=`dim 3`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(mxylimp3,xylimp3))  + scale_y_continuous(limits=c(mxylimp3,xylimp3))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 2 (",MFA2, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none");
	p3 <- p3+ geom_segment(aes(x=x , y=x, xend=`dim 2`, yend=`dim 3`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=Df5plot$colPalette,labels=Df5plot$labels) +ggforce::geom_circle(inherit.aes = FALSE, aes(x0=0, y0=0,r=xylimp3-0.4))+coord_fixed()
	#plot" label
	xylimp3Label=max(c(Df5plot[,2],Df5plot[,3]))+0.7 + margin
	mxylimp3Label=min(c(Df5plot[,2],Df5plot[,3]))-0.7 - margin

	p3label <- ggplot(Df5plot, aes(x=`dim 2`, y=`dim 3`,color=labels)) + geom_point(color="#f0f0f5")  + 
	scale_x_continuous(limits=c(mxylimp3Label,xylimp3Label))  + scale_y_continuous(limits=c(mxylimp3Label,xylimp3Label))+ 
	theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ labs(x = paste0("dim 2 (",MFA2, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none");
	p3label <- p3label+ geom_segment(aes(x=x , y=x, xend=`dim 2`, yend=`dim 3`),
	                  arrow=arrow(length=unit(0.3,"cm")))+		 	  
	     scale_color_manual(values=Df5plot$colPalette,labels=Df5plot$labels) +ggforce::geom_circle(inherit.aes = FALSE, aes(x0=0, y0=0,r=xylimp3Label-0.15))+coord_fixed()
	p3labelShort <- p3label + geom_text(aes(x=`dim 2`, y=`dim 3`,label=labelsParameter,size=0.1,vjust = 0.7,hjust=-0.13, angle=45))
	p3labelFull <- p3label + geom_text(aes(x=`dim 2`, y=`dim 3`,label=labels,size=0.1,vjust = 0.7,hjust=-0.13, angle=45))
	p3label_nb <- p3label + geom_text(aes(x=`dim 2`, y=`dim 3`,label=nb,size=0.1,vjust = ifelse(`dim 3`>0, -.2,`dim 3`+2.4),hjust=ifelse(`dim 3`>0, `dim 2`-0.1,`dim 2`+0.8), angle=45)) #ifelse(`dim 1`>0, 0.3,`dim 1`-0.5)

	#legend
	pc_qp6<-qplot(data=Df5plot, x=`dim 2`, y=`dim 3`,color=labels )+
		 	 geom_point(size=4,stat = "identity") + geom_hline(yintercept = 0, colour = "gray65") +
		 	  geom_vline(xintercept = 0, colour = "gray65")+ 
		 	  scale_color_manual(values=Df5plot$colPalette,labels=Df5plot$labels)+
		 	    labs(x = paste0("dim 2 (",MFA2, "% )"), 
		 	    	y=paste0("dim 3 (",MFA3, "% )")) +
		 	    guides(colour = guide_legend("labels",nrow = nrowLegend, byrow = TRUE))+ 
		 	    theme_gray(base_size = 12)+
		 	    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
		 	legend.title = element_text(size=12, face="bold"),legend.position="bottom", 
		 	legend.box = "horizontal",legend.text=element_text(size=10,face="bold"),
		 	legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line"));

	pc_qp6_nb<-qplot(data=Df5plot, x=`dim 2`, y=`dim 3`,color=labels )+
		 	 geom_point(size=4,stat = "identity") + geom_hline(yintercept = 0, colour = "gray65") +
		 	  geom_vline(xintercept = 0, colour = "gray65")+ 
		 	  scale_color_manual(values=Df5plot$colPalette,labels=Df5plot$labels_nb)+
		 	    labs(x = paste0("dim 2 (",MFA2, "% )"), 
		 	    	y=paste0("dim 3 (",MFA3, "% )")) +
		 	    guides(colour = guide_legend("labels",nrow = nrowLegend, byrow = TRUE))+ 
		 	    theme_gray(base_size = 12)+
		 	    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
		 	legend.title = element_text(size=12, face="bold"),legend.position="bottom", 
		 	legend.box = "horizontal",legend.text=element_text(size=10,face="bold"),
		 	legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line"));

	#labels
	upText <-paste0("    ",nameExperiment,": input file with ",nbVarTotalF," phenotypic variables");
	upText.p <- ggparagraph(text = upText, face = "bold", size = 14, color = "black");
	midText <-paste0("    ",nameExperiment,": MFA considering ",nbVarF," non-correlated variables");
	midText.p <- ggparagraph(text = midText, face = "bold", size = 14, color = "black");

	# Arrange the plots on the same page  --> rats and all other plots x11(height=5, width=10.2)
	#X11(height=11.3, width=8.2) #8.3 x 11.7 i
	legend <- cowplot::get_legend(pc_qp6); 
	legend_nb <- cowplot::get_legend(pc_qp6_nb); 
	blank<-rectGrob(gp=gpar(col="white")) # make a white spacer grob
	p1a <- ggarrange(upText.p, ncol = 1, nrow = 1);
	p1b <- ggarrange(midText.p, ncol = 1, nrow = 1);
	p2a <- ggarrange(p1, p2, p3,  ncol = 3, nrow = 1, heights = c(1, 1, 1));
	p12 <- ggarrange(p1a,p1b, p2a, ncol = 1, nrow = 3,heights = c(0.1,0.1, 1));
	p4a <- ggarrange(p1labelShort, p2labelShort, p3labelShort +  theme(legend.position="none"),  ncol = 3, nrow = 1, heights = c(1, 1, 1));
	p4b <- ggarrange(p1labelFull, p2labelFull, p3labelFull +  theme(legend.position="none"),  ncol = 3, nrow = 1, heights = c(1, 1, 1));
	p4_nb <- ggarrange(p1label_nb, p2label_nb, p3label_nb +  theme(legend.position="none"),  ncol = 3, nrow = 1, heights = c(1, 1, 1));

    p2v2a <-ggarrange(p1labelShort, p2labelShort, p3labelShort +  theme(legend.position="none"),  ncol = 1, nrow = 3, widths = c(1, 1, 1));
    p2v2b <-ggarrange(p1labelFull, p2labelFull, p3labelFull +  theme(legend.position="none"),  ncol = 1, nrow = 3, widths = c(1, 1, 1));
    p2v2_nb <-ggarrange(p1label_nb, p2label_nb, p3label_nb +  theme(legend.position="none"),  ncol = 1, nrow = 3, widths = c(1, 1, 1));
    p2v2c <-ggarrange(p1, p2, p3 +  theme(legend.position="none"),  ncol = 1, nrow = 3, widths = c(1, 1, 1));
	p12v2a <- ggarrange(p1a,p1b,p2v2a,legend, ncol = 1, nrow = 4,heights = c(0.02,0.02,1,0.5));
	p12v2b <- ggarrange(p1a,p1b,p2v2b,legend, ncol = 1, nrow = 4,heights = c(0.02,0.02,1,0.5));
	p12v2c <- ggarrange(p1a,p1b,p2v2c,legend, ncol = 1, nrow = 4,heights = c(0.02,0.02,1,0.5));
	p12v2_nb <- ggarrange(p1a,p1b,p2v2_nb,legend_nb, ncol = 1, nrow = 4,heights = c(0.02,0.02,1,0.5));

	#individual labelled plots
	p1labelarrangedA <- ggarrange(blank,p1a,p1b,p1labelShort, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p2labelarrangedA <- ggarrange(blank,p1a,p1b,p2labelShort, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p3labelarrangedA <- ggarrange(blank,p1a,p1b,p3labelShort, ncol = 1, nrow = 4, heights = c(0.01,0.05,0.05,1));
	p1labelarrangedB <- ggarrange(blank,p1a,p1b,p1labelFull, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p2labelarrangedB <- ggarrange(blank,p1a,p1b,p2labelFull, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p3labelarrangedB <- ggarrange(blank,p1a,p1b,p3labelFull, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p1labelarranged_nb <- ggarrange(blank,p1a,p1b,p1label_nb, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p2labelarranged_nb <- ggarrange(blank,p1a,p1b,p2label_nb, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	p3labelarranged_nb <- ggarrange(blank,p1a,p1b,p3label_nb, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));

	#final grouped objects
	pf1a <- ggarrange(blank,p12,blank,p4a,legend, ncol = 1,nrow = 5,heights = c(0.1,1,0.1,1,0.5)); #with/wo mabels
	pf1b <- ggarrange(blank,p12,blank,p4b,legend, ncol = 1,nrow = 5,heights = c(0.1,1,0.1,1,0.5)); #with/wo mabels
	pf1_nb <- ggarrange(blank,p12,blank,p4_nb,legend_nb, ncol = 1,nrow = 5,heights = c(0.1,1,0.1,1,0.5)); #with/wo mabels
	pf2 <- ggarrange(blank,p12,legend, ncol = 1,nrow = 3,heights = c(0.1,1,0.5)); #just wo labels
	p4v2a <-	ggarrange(blank,p12v2a, ncol = 1,nrow = 2,heights = c(0.05,1)); #just with labels
	p4v2b <-	ggarrange(blank,p12v2b, ncol = 1,nrow = 2,heights = c(0.05,1)); #just with labels
	p4v2c <-	ggarrange(blank,p12v2c, ncol = 1,nrow = 2,heights = c(0.05,1)); #just with labels
	p4v2_nb <-	ggarrange(blank,p12v2_nb, ncol = 1,nrow = 2,heights = c(0.05,1)); #just with labels

	#saving into pdf final plots arranged
	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot5,nameExperiment, "_quantitativeVarCoordinates_IndividualPlot.pdf"),height=7, width=11)
	print(p1labelarrangedA);
	print(p2labelarrangedA);
	print(p3labelarrangedA);
	print(p1labelarrangedB);
	print(p2labelarrangedB);
	print(p3labelarrangedB);
	print(p1labelarranged_nb);
	print(p2labelarranged_nb);
	print(p3labelarranged_nb);
	dev.off();

	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot5,nameExperiment, "_quantitativeVarCoordinates_labelGrouped",".pdf"),height=20, width=12)
	print(p4v2a);
	print(p4v2b);
	print(p4v2c);
	print(p4v2_nb);
	dev.off();
	
	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot5,nameExperiment, "_quantitativeVarCoordinates_unlabelGrouped",".pdf"),height=8, width=14)
	print(pf2);
	dev.off();

	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot5,nameExperiment, "_quantitativeVarCoordinates_allGrouped",".pdf"),height=18, width=20)
	print(pf1a);
	print(pf1b);
	print(pf1_nb);
	dev.off();
};




#' plot6 -	Quantitative variables distribution in PCA
#'
#' @name plot6_qualitative_vars_discrimination_PCA.function 
#'
#' @param name name for the plot or specific to the session
#'
#' @param catVarName Column name of the qualitative
#'  variable or categorical variable that is your dependant 
#'  variable, the one for which one you want to identify the most discriminative variables  
#' 
#' @param nameExperiment The name of the experiment together with the number of variables we included. 
#' if is all model, allthe variables. Instead if is the selected model is less variables. 
#' ex: 
#' 1. nameExperiment <- paste0(name,"_", length(colnames(phenoSocial_4analysis)), " variables"); #all variables
#' 2. nameExperiment <- paste0(name,"_", length(colnames(phenoAll_4analysis_woHighCor_sel1)), " variables"); · without highly correlated
#' 3. nameExperiment <- paste0(name,"_", " more informative variables ",length(colnames(parameters_more_30perc_informative_Genotype)) ) #selected variables
#'
#' @param nrowLegend Number of rows the legend should be split to, to format properly for publications.
#'
#' @param Cat_variableToShow  name of the specific categorical variable to show. Specify isntead  "All"
#' to show all in the plot.
#'
#' @param nameModel Short naming for the model you are using.
#' Ex: "fullMdl" for the model including all variables;"Sel1" for the model not including the
#' highly correlated variables, "Sel30" for the model including only the variables contributing
#' more thana 30%
#' 
#' @return The plot 6 is generated in the plot6 folder in three different kind of organizations
#' to select the better for publications (landscape, vertical, 1 plot alone or all plots together.
#'
#' @examples
#' # show the "Genotype" variable
#' plot6_qualitative_vars_discrimination_PCA.function(name,catVarName,nameExperiment,4,"Genotype",nameModel);
#' # show all the categorical variables . Ex. Sex and Genotype
#' plot6_qualitative_vars_discrimination_PCA.function(name,catVarName,nameExperiment,4,"All",nameModel);
#'
#' @export
plot6_qualitative_vars_discrimination_PCA.function <- function(name,catVarName,nameExperiment,nrowLegend,Cat_variableToShow,nameModel){
    load( file=paste0(wd, f_pheno_discrimination,f_results, f_Rdata,name,"_",nameModel,"_",catVarName,"_Robject_for_plots.RData")) #namesVariablesDf, df1, df2, df3,df4, df5, df6, MFA1, MFA2, MFA3,nbVarF, nbVarTotalF,
	Df6plot <- df6

	if (Cat_variableToShow=="All"){
		Df6plot <- Df6plot;
	} else {
		Df6plot <- Df6plot[grep(Cat_variableToShow,Df6plot$Variable),];
	};

	Df6plot$labels <- factor(Df6plot$labels, levels=unique(Df6plot$labels))
	
	#x and y axis limits
	min_value_3d <- min(c(Df6plot[,1],Df6plot[,2],Df6plot[,3])); 
	max_value_3d <- max(c(Df6plot[,1],Df6plot[,2],Df6plot[,3])); 

	#2) plots for the variables selected on the analysis
	pc_qp1<- qplot(x=`dim 1`, y=`dim 2`, data=as.data.frame(Df6plot), color=labels) + scale_x_continuous(limits=c(min_value_3d -0.05,max_value_3d +0.02))  + scale_y_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.05)) + theme_gray(base_size = 11) + theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ scale_color_manual(breaks = labels, values=Df6plot$colPalette)+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 2 (",MFA2, "% )"))+ theme(legend.position="none");  #take out the ,labels=my_leg1 to not put the italics
	pc_qp1label<-pc_qp1 +geom_text(aes(label=labels,size=2,vjust = -0.75,hjust=0.4,angle=10))  #vjust = 0.1,

	pc_qp2<-qplot(x=`dim 1`, y=`dim 3`, data=as.data.frame(Df6plot), color=labels) + scale_x_continuous(limits=c(min_value_3d -0.05,max_value_3d +0.02))  + scale_y_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.05)) + theme_gray(base_size = 11)+ theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ scale_color_manual(breaks = labels, values=Df6plot$colPalette)+ labs(x = paste0("dim 1 (",MFA1, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none"); 
	pc_qp2label<-pc_qp2 +geom_text(aes(label=labels,size=2,vjust = -0.75,hjust=0.4,angle=10))  #vjust = 0.1,

	pc_qp3<-qplot(x=`dim 2`, y=`dim 3`, data=as.data.frame(Df6plot), color=labels) + scale_x_continuous(limits=c(min_value_3d -0.05,max_value_3d +0.02))  + scale_y_continuous(limits=c(min_value_3d -0.02,max_value_3d +0.05)) + theme_gray(base_size = 11)+ theme(axis.text=element_text(size=14), axis.title=element_text(size=15,face="bold"),legend.title = element_text(size=9, face="bold"),legend.position="bottom", legend.box = "horizontal",legend.text=element_text(size=9),legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line")) + geom_point(size=4) + geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ scale_color_manual(breaks = labels, values=Df6plot$colPalette)+  labs(x = paste0("dim 2 (",MFA2, "% )"), y=paste0("dim 3 (",MFA3, "% )"))+ theme(legend.position="none"); 
	pc_qp3label<-pc_qp3 +geom_text(aes(label=labels,size=2,vjust = -0.75,hjust=0.4,angle=10))  #vjust = 0.1,

	pc_qp6<-qplot(data=Df6plot, x=`dim 2`, y=`dim 3`,color=labels )+
		 	 geom_point(size=4,stat = "identity") + geom_hline(yintercept = 0, colour = "gray65") +
		 	  geom_vline(xintercept = 0, colour = "gray65")+ 
		 	  scale_color_manual(values=Df6plot$colPalette,labels=Df6plot$labels)+
		 	    labs(x = paste0("dim 2 (",MFA2, "% )"), 
		 	    	y=paste0("dim 3 (",MFA3, "% )")) +
		 	    guides(colour = guide_legend("labels",nrow = nrowLegend, byrow = TRUE))+ 
		 	    theme_gray(base_size = 12)+
		 	    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
		 	legend.title = element_text(size=12, face="bold"),legend.position="bottom", 
		 	legend.box = "horizontal",legend.text=element_text(size=10,face="bold"),
		 	legend.key.width=unit(1.1,"line"),legend.key.height=unit(1.1,"line"));

	#labels
	upText <-paste0("    ",nameExperiment,": input file with ",nbVarTotalF," phenotypic variables");
	upText.p <- ggparagraph(text = upText, face = "bold", size = 14, color = "black");
	midText <-paste0("    ",nameExperiment,": MFA considering ",nbVarF," non-correlated variables");
	midText.p <- ggparagraph(text = midText, face = "bold", size = 14, color = "black");


	# Arrange the plots on the same page  --> rats and all other plots x11(height=5, width=10.2)
	#X11(height=13.3, width=20.2) #8.3 x 11.7 i
	legend <- cowplot::get_legend(pc_qp6); 
	blank<-rectGrob(gp=gpar(col="white")) # make a white spacer grob
	p1a <- ggarrange(upText.p, ncol = 1, nrow = 1);
	p1b <- ggarrange(midText.p, ncol = 1, nrow = 1);
	p2 <- ggarrange(pc_qp1, pc_qp2, pc_qp3,  ncol = 3, nrow = 1, heights = c(1, 1, 1));
	p12 <- ggarrange(p1a,p1b, p2, ncol = 1, nrow = 3,heights = c(0.1,0.1, 1));
	p4 <- ggarrange(pc_qp1label, pc_qp2label, pc_qp3label +  theme(legend.position="none"),  ncol = 3, nrow = 1, heights = c(1, 1, 1));

    p2v2 <-ggarrange(pc_qp1label, pc_qp2label, pc_qp3label +  theme(legend.position="none"),  ncol = 1, nrow = 3, widths = c(1, 1, 1));
	p12v2 <- ggarrange(p1a,p1b,p2v2, ncol = 1, nrow = 3,heights = c(0.02,0.02,1));
	#p4v2 <-	ggarrange(blank,p12v2,legend, ncol = 1,nrow = 3,heights = c(0.1,1,0.5));w legend
	#p4v2 <-	ggarrange(blank,p12v2, ncol = 1,nrow = 2,heights = c(0.05,1)); #wo legend

	#individual labelled plots
	pc_qp1labelarranged <- ggarrange(blank,p1a,p1b,pc_qp1label, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	pc_qp2labelarranged <- ggarrange(blank,p1a,p1b,pc_qp2label, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));
	pc_qp3labelarranged <- ggarrange(blank,p1a,p1b,pc_qp3label, ncol = 1, nrow = 4,heights = c(0.01,0.05,0.05,1));

	#final grouped objects
	pf1 <- ggarrange(blank,p12,blank,p4,legend, ncol = 1,nrow = 5,heights = c(0.1,1,0.1,1,0.5)); #with/wo mabels
	pf2 <- ggarrange(blank,p12,legend, ncol = 1,nrow = 3,heights = c(0.1,1,0.5)); #just wo labels
	p4v2 <-	ggarrange(blank,p12v2, ncol = 1,nrow = 2,heights = c(0.05,1)); #just with labels

	# saving to pdf
	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot6,nameExperiment,"_",Cat_variableToShow, "_levelsComponents_sex_Geno_labelGrouped",".pdf"),height=12, width=7)
	print(p4v2);
	dev.off();
	
	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot6,nameExperiment,"_",Cat_variableToShow, "_levelsComponents_sex_Geno_unlabelGrouped",".pdf"),height=6, width=12)
	print(pf2);
	dev.off();

	pdf(paste0(wd,f_pheno_discrimination, f_results,f_plot6,nameExperiment,"_",Cat_variableToShow, "_levelsComponents_sex_Geno_allGrouped",".pdf"),height=11, width=14)
	print(pf1);
	dev.off();
};




scatterplot3d.function <- function(input_scatterplot3d,name, catVarName, nameExperiment,colorByColumn,colPalette_toPlot,height,width,m,ncolLegend,nameModel){
	#library("scatterplot3d");
    load( file=paste0(wd, f_pheno_discrimination,f_results, f_Rdata,name,"_",nameModel,"_",catVarName,"_Robject_for_plots.RData")) 
    #namesVariablesDf, df1, df2, df3,df4, df5, df6, MFA1, MFA2, MFA3,nbVarF, nbVarTotalF,

    #pre-processing of the input file to reorder cat variables first
	quantVariables <- input_scatterplot3d[,sapply(input_scatterplot3d, is.numeric)];
	colVariablesNames <- colnames(input_scatterplot3d[,grep(colPalette_toPlot,colnames(input_scatterplot3d))]);
	quantVarsNames <- colnames(quantVariables);
	qualiVariablesNames<- colnames(input_scatterplot3d[ , !(colnames(input_scatterplot3d) %in% quantVarsNames)]);
	qualiVarsNumber<- length(qualiVariablesNames);
	catVarsNumber<- length(quantVariables);

	input_scatterplot3dF<- input_scatterplot3d[ , c(qualiVariablesNames,quantVarsNames)];

	for (i in 1:qualiVarsNumber){
		input_scatterplot3dF[,i]<-factor(input_scatterplot3dF[,i],levels=unique(input_scatterplot3dF[,i])); 
		i <- i+1;
	};

	input_scatterplot3dF$colorByColumn <- gsub("_"," ",input_scatterplot3dF[,colorByColumn]);
	input_scatterplot3dF$colorByColumn<-factor(input_scatterplot3dF$colorByColumn,levels=unique(input_scatterplot3dF$colorByColumn)); 

	input_scatterplot3dF$colPalette_toPlot <- input_scatterplot3dF[,colPalette_toPlot];
	#input_scatterplot3dF$colPalette_toPlot<-factor(input_scatterplot3dF$colPalette_toPlot,levels=unique(input_scatterplot3dF$colPalette_toPlot)); 
	##################################################################################################
	plot4legend <- qplot(data=input_scatterplot3dF, x=`dim 2`, y=`dim 3`,color=colPalette_toPlot )+
		geom_point(size=4,stat = "identity") + geom_hline(yintercept = 0, colour = "gray65") +
		geom_vline(xintercept = 0, colour = "gray65")+ 
		scale_color_manual(values=levels(input_scatterplot3dF$colPalette_toPlot),labels=unique(input_scatterplot3dF$colorByColumn))+
		guides(colour = guide_legend("labels",ncol = ncolLegend, byrow = FALSE))+ 
		theme_gray(base_size = 12)+ 
		theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"),
		 	legend.title = element_blank(),legend.position="bottom", 
		 	legend.box = "horizontal",legend.text=element_text(size=14),
		 	legend.key.width=unit(1.4,"line"),legend.key.height=unit(1.4,"line")
		);

	legend <- cowplot::get_legend(plot4legend); 

	x11(height=height, width=width);
	layout(m)

	par(mar=c(2.5, 5.5, 2, 4));
	maxEGs  <-  max(input_scatterplot3dF[,c("dim 1", "dim 2", "dim 3")]) + 0.2;
	minEGs  <-  min(input_scatterplot3dF[,c("dim 1", "dim 2", "dim 3")]) - 0.2;

	scatterplot3d(input_scatterplot3dF$`dim 1`,input_scatterplot3dF$`dim 2`,input_scatterplot3dF$`dim 3`,       
	                color=input_scatterplot3dF$colPalette_toPlot, pch=16, 
	                type="p",cex.symbols=3.1,angle=35, col.axis="grey28", bg="grey81",lty.hide=2, lab = c(2,2), lab.z = 2,
	                cex.axis=1.1,cex.lab=0.9, y.margin.add=0.4, xlim=c(minEGs,maxEGs), ylim=c(minEGs,maxEGs), zlim=c(minEGs,maxEGs),
	                col.grid="lightblue", 
	                main=paste0(nameExperiment," \n (",nbVarF, " selected variables from ", nbVarTotalF," )"),
	                xlab=paste0("Dim 1 (",MFA1,"%)"),
	                ylab=paste0("                           Dim 2 (",MFA2,"%)"),
	                zlab=paste0("Dim 3 (",MFA3,"%)"))

	frame();
	vps <- baseViewports();
	pushViewport(vps$inner, vps$figure, vps$plot);
	#
	grid.draw(legend);
	dev.copy2pdf(device = quartz, file = paste0(wd,f_pheno_discrimination, f_results,f_plot3,gsub(" ","_",nameExperiment),"_",
		catVarName,"_",colorByColumn,"_",colPalette_toPlot, "_3dpca.pdf"),onefile=TRUE);

	dev.off();

};


#' parallel plot
#'
#' @name paralellePlot.function 
#'
#' @param inputDf_parallele_Sumup Input files for computing the parallel plot name using 
#' the mean and sd and CI of the observations. 
#'
#' @param inputDf_parallele_Ind Input files for computing the parallel plot per observation
#' ' Column name of the qualitative
#' 
#' @param parallelePlot_phenoColourDf Dataframe containing at least 2 columns, 
#' 1: "condition", containing the possible categories. Ex Genotype het, wt, homo
#' 2: "colPalette_cond_Geno_Sex" hex colour associated to each variable. 
#' Ex: for condition="wt_female ", use colPalette="#800080" 
#' 
#' @param nameOutput Suffix added to the name of the plot, specfic to this project 
#'
#' @param height  Height in inches to save the plot
#'
#' @param width Width in inches to save the plot
#' Ex: "fullMdl" for the model including all variables;"Sel1" for the model not including the
#' highly correlated variables, "Sel30" for the model including only the variables contributing
#' more thana 30%
#' 
#' @return The plot 6 is generated in the plot6 folder in three different kind of organizations
#' to select the better for publications (landscape, vertical, 1 plot alone or all plots together.
#'
#' @examples
#height <-4
#width <-4.5
#paralellePlot.function(inputDf_parallele_Sumup_Genotype,inputDf_parallele_Ind,parallelePlot_phenoColourDf,"AllGenotypes_Sexes",height,width)

#' @export
paralellePlot.function <- function(inputDf_parallele_Sumup,inputDf_parallele_Ind,parallelePlot_phenoColourDf,nameOutput,height,width){
	##########################################@
	#Two options to built the paralelle plot:
	##########################################@
	#1.paralel plot with the mean + stats ( sd or error or CI) using inputDf_parallele
	#2. paralel plot with all the individual values -> this is too messy and not informative thus not implementing it using inputDf_parallele
	##########################################@
	inputDf_parallele_mean <- inputDf_parallele_Sumup[,c(grep("condition",colnames(inputDf_parallele_Sumup)),grep("mean", colnames(inputDf_parallele_Sumup))) ];
	inputDf_parallele4Plot_plot_mean<-inputDf_parallele_mean;
	inputDf_parallele4Plot_plot_mean <- inputDf_parallele4Plot_plot_mean[ , colSums(is.na(inputDf_parallele4Plot_plot_mean)) == 0];
	inputDf_parallele4Plot_plot_mean<- left_join(inputDf_parallele4Plot_plot_mean, parallelePlot_phenoColourDf,by="condition");
	nColPlot_mean <-2:(length(unique(colnames(inputDf_parallele4Plot_plot_mean)))-1);
	inputDf_parallele4Plot_plot_mean$colPalette <- inputDf_parallele4Plot_plot_mean[,grep("colPalette",colnames(inputDf_parallele4Plot_plot_mean))];
	inputDf_parallele4Plot_plot_mean$condition <- factor(inputDf_parallele4Plot_plot_mean$condition, levels=unique(inputDf_parallele4Plot_plot_mean$condition));
	###################################################################@ 
	#1) parallel plot using the means per parameter
	###################################################################@

	plot_noScaling.tendency <- ggparcoord(inputDf_parallele4Plot_plot_mean,
	    columns = nColPlot_mean, groupColumn = 1, order = "allClass",
	    scale="globalminmax",
	    showPoints = FALSE, 
	    title = "No scTRUE",
	    alphaLines = 1) + coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
		scale_color_viridis(discrete=TRUE) + theme_ipsum()+scale_alpha(range = c(0.3, 1))+
		theme(legend.position="bottom",plot.title = element_text(size=13),legend.key.size = unit(1,"line"),legend.text=element_text(size=14)) + xlab("") +
     	guides(color = guide_legend(override.aes = list(size = 1.7)),fill=guide_legend(nrow=3, byrow = TRUE) )+ scale_color_manual(values=inputDf_parallele4Plot_plot_mean$colPalette);
	
	plot_noScaling_points.tendency <- ggparcoord(inputDf_parallele4Plot_plot_mean,
	    columns = nColPlot_mean, groupColumn = 1, order = "allClass",
	    scale="globalminmax",
	    showPoints = TRUE, 
	    title = "No scTRUE",
	    alphaLines = 0) + coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
		scale_color_viridis(discrete=TRUE) + theme_ipsum()+scale_alpha(range = c(0.3, 1))+
		theme(legend.position="bottom",plot.title = element_text(size=13),legend.key.size = unit(1,"line"),legend.text=element_text(size=14)) + xlab("") +
     	guides(color = guide_legend(override.aes = list(size = 1.7)),fill=guide_legend(nrow=3, byrow = TRUE) )+ scale_color_manual(values=inputDf_parallele4Plot_plot_mean$colPalette);


	plot_Standard.tendency<- ggparcoord(inputDf_parallele4Plot_plot_mean,
	    columns = nColPlot_mean, groupColumn = 1, order = "allClass",
	    scale="uniminmax",
	    showPoints = FALSE, 
	    title = "Standardize to Min = 0 and Max = 1",
	    alphaLines = 1) + coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
		scale_color_viridis(discrete=TRUE) + theme_ipsum()+scale_alpha(range = c(0.3, 1))+
		theme(legend.position="bottom",plot.title = element_text(size=13),legend.key.size = unit(1,"line"),legend.text=element_text(size=14)) + xlab("") +
     	guides(color = guide_legend(override.aes = list(size = 1.7)),fill=guide_legend(nrow=3, byrow = TRUE) )+ scale_color_manual(values=inputDf_parallele4Plot_plot_mean$colPalette);


	plot_StandardCenter.tendency<-ggparcoord(inputDf_parallele4Plot_plot_mean,
	    columns = nColPlot_mean, groupColumn = 1, order = "allClass",
	    scale="center",
	    showPoints = FALSE, 
	    title = "Standardize and center variables",
	    alphaLines = 1) + coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
		scale_color_viridis(discrete=TRUE) + theme_ipsum()+scale_alpha(range = c(0.3, 1))+
		theme(legend.position="bottom",plot.title = element_text(size=13),legend.key.size = unit(1,"line"),legend.text=element_text(size=14)) + xlab("") +
     	guides(color = guide_legend(override.aes = list(size = 1.7)),fill=guide_legend(nrow=3, byrow = TRUE) )+ scale_color_manual(values=inputDf_parallele4Plot_plot_mean$colPalette);

	plot_Univariately.tendency<-ggparcoord(inputDf_parallele4Plot_plot_mean,
	    columns = nColPlot_mean, groupColumn = 1, order = "allClass",
	    scale="std",
	    showPoints = FALSE, 
	    title = "Normalize univariately (substract mean & divide by sd)",
	    alphaLines = 1) + coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
		scale_color_viridis(discrete=TRUE) + theme_ipsum()+scale_alpha(range = c(0.3, 1))+
		theme(legend.position="bottom",plot.title = element_text(size=13),legend.key.size = unit(1,"line"),legend.text=element_text(size=14)) + xlab("") +
     	guides(color = guide_legend(override.aes = list(size = 1.7)),fill=guide_legend(nrow=3, byrow = TRUE) )+ scale_color_manual(values=inputDf_parallele4Plot_plot_mean$colPalette);

	###################################################################@ 
	#2) parallele plot showing all values per parameter
	###################################################################@
	inputDf_parallele_Ind<- left_join(inputDf_parallele_Ind, parallelePlot_phenoColourDf,by="condition");
	nColPlot_Ind <-2:(length(unique(colnames(inputDf_parallele_Ind)))-1);

	plot_noScaling.ind <- ggparcoord(inputDf_parallele_Ind,
	    columns = nColPlot_Ind, groupColumn = 1, order = "allClass",
	    scale="globalminmax",
	    showPoints = TRUE, 
	    title = "No scTRUE all Ind",
	    alphaLines = 0, mapping=aes(color=as.factor(condition))) + coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+ 
		scale_color_viridis(discrete=TRUE) + theme_ipsum()+
		theme(legend.position="bottom",plot.title = element_text(size=13),legend.key.size = unit(1,"line"),legend.text=element_text(size=14)) + xlab("") +
     	guides(color = guide_legend(override.aes = list(size = 1.7)),fill=guide_legend(nrow=2, byrow = TRUE) )+ scale_color_manual(values=unique(inputDf_parallele_Ind$colPalette));
	

	plot_Standard.ind<- ggparcoord(inputDf_parallele_Ind,
	    columns = nColPlot_Ind, groupColumn = 1, order = "allClass",
	    scale="uniminmax",
	    showPoints = TRUE, 
	    title = "Standardize to Min = 0 and Max = 1 all Ind",
	    alphaLines = 0, mapping=aes(color=as.factor(condition))) + coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
		scale_color_viridis(discrete=TRUE) + theme_ipsum()+
		theme(legend.position="bottom",plot.title = element_text(size=13),legend.key.size = unit(1,"line"),legend.text=element_text(size=14)) + xlab("") +
     	guides(color = guide_legend(override.aes = list(size = 1.7)),fill=guide_legend(nrow=2, byrow = TRUE) )+ scale_color_manual(values=unique(inputDf_parallele_Ind$colPalette));

	plot_StandardCenter.ind<-ggparcoord(inputDf_parallele_Ind,
	    columns = nColPlot_Ind, groupColumn = 1, order = "allClass",
	    scale="center",
	    showPoints = TRUE, 
	    title = "Standardize and center variables all Ind",
	    alphaLines = 0, mapping=aes(color=as.factor(condition))) + coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
		scale_color_viridis(discrete=TRUE) + theme_ipsum()+
		theme(legend.position="bottom",plot.title = element_text(size=13),legend.key.size = unit(1,"line"),legend.text=element_text(size=14)) + xlab("") +
     	guides(color = guide_legend(override.aes = list(size = 1.7)),fill=guide_legend(nrow=2, byrow = TRUE) )+ scale_color_manual(values=unique(inputDf_parallele_Ind$colPalette));


	###################################################################@ 
	#3) adding CI in the plot using a line up and down in softer color
	###################################################################@
	inputDf_parallele_mean_plot2 <- inputDf_parallele_Sumup[,c(grep("condition",colnames(inputDf_parallele_Sumup)),grep("mean", colnames(inputDf_parallele_Sumup))) ];
	inputDf_parallele_ci_lower <- inputDf_parallele_Sumup[,c(grep("condition",colnames(inputDf_parallele_Sumup)),grep("ci_lower", colnames(inputDf_parallele_Sumup))) ];
	inputDf_parallele_ci_upper <- inputDf_parallele_Sumup[,c(grep("condition",colnames(inputDf_parallele_Sumup)),grep("ci_upper", colnames(inputDf_parallele_Sumup))) ];

	colnames(inputDf_parallele_mean_plot2) <- gsub("mean ","",colnames(inputDf_parallele_mean_plot2));
	colnames(inputDf_parallele_ci_lower) <- gsub("ci_lower ","",colnames(inputDf_parallele_ci_lower));
	colnames(inputDf_parallele_ci_upper) <- gsub("ci_upper ","",colnames(inputDf_parallele_ci_upper));

	inputDf_parallele_mean_plot2$alpha <- 1
	inputDf_parallele_ci <- bind_rows(inputDf_parallele_ci_lower,inputDf_parallele_ci_upper);
	inputDf_parallele_ci$alpha <- 0.3
	inputDf_parallele_ci <- bind_rows(inputDf_parallele_mean_plot2,inputDf_parallele_ci);
	nColPlot <- 2:(length(unique(colnames(inputDf_parallele_ci)))-3);
	inputDf_parallele_ci<- left_join(inputDf_parallele_ci, parallelePlot_phenoColourDf,by="condition");
	inputDf_parallele_ci$colPalette <- inputDf_parallele_ci[,grep("colPalette",colnames(inputDf_parallele_ci))]
	#inputDf_parallele_ci[,"colPalette"] <- factor(inputDf_parallele_ci[,"colPalette"],levels=unique(inputDf_parallele_ci[,"colPalette"]));
	inputDf_parallele_ci$condition <- factor(inputDf_parallele_ci$condition, levels=unique(inputDf_parallele_ci$condition));
	                                                     
	#erorr here!
	plot_noScaling.tendency.CI <- ggparcoord(inputDf_parallele_ci, columns = nColPlot, groupColumn = 1, order = "allClass",alpha ="alpha", scale="globalminmax",showPoints = FALSE, title = "No scTRUE + CI") + 
			coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
		scale_color_viridis(discrete=TRUE) + theme_ipsum()+scale_alpha(range = c(0.3, 1))+
		theme(legend.position="bottom",plot.title = element_text(size=13),legend.key.size = unit(1,"line"),legend.text=element_text(size=14)) + xlab("") +
     	guides(color = guide_legend(override.aes = list(size = 1.7)),fill=guide_legend(nrow=2, byrow = TRUE) )+ scale_color_manual(values=inputDf_parallele_ci$colPalette);


	plot_Standard.tendency.CI<- ggparcoord(inputDf_parallele_ci,
	    columns = nColPlot, groupColumn = 1, order = "allClass",alpha ="alpha",
	    scale="uniminmax", showPoints = FALSE, title = "Standardize to Min = 0 and Max = 1 + CI") + 
		coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
		scale_color_viridis(discrete=TRUE) + theme_ipsum()+scale_alpha(range = c(0.3, 1))+
		theme(legend.position="bottom",plot.title = element_text(size=13),legend.key.size = unit(1,"line"),legend.text=element_text(size=14)) + xlab("") +
     	guides(color = guide_legend(override.aes = list(size = 1.7)),fill=guide_legend(nrow=2, byrow = TRUE) )+ scale_color_manual(values=inputDf_parallele_ci$colPalette);


	plot_StandardCenter.tendency.CI<-ggparcoord(inputDf_parallele_ci,
	    columns = nColPlot, groupColumn = 1, order = "allClass",alpha ="alpha",
	    scale="center", showPoints = FALSE, title = "Standardize and center variables + CI") + 
		coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
	    scale_color_viridis(discrete=TRUE) + theme_ipsum()+
		theme(legend.position="bottom",plot.title = element_text(size=13),legend.key.size = unit(1,"line"),legend.text=element_text(size=14)) + xlab("") +
     	guides(color = guide_legend(override.aes = list(size = 1.7)),fill=guide_legend(nrow=2, byrow = TRUE) )+ scale_color_manual(values=inputDf_parallele_ci$colPalette);

	plot_Univariately.tendency.CI<-ggparcoord(inputDf_parallele_ci,
	    columns = nColPlot, groupColumn = 1, order = "allClass",alpha ="alpha",
	    scale="std", showPoints = FALSE, title = "Normalize univariately (substract mean & divide by sd) + CI") +
		coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
	    scale_color_viridis(discrete=TRUE) + theme_ipsum()+
		theme(legend.position="bottom",plot.title = element_text(size=13),legend.key.size = unit(1,"line"),legend.text=element_text(size=14)) + xlab("") +
     	guides(color = guide_legend(override.aes = list(size = 1.7)),fill=guide_legend(nrow=2, byrow = TRUE) )+ scale_color_manual(values=inputDf_parallele_ci$colPalette);


    #save the plots

	#extrafont::loadfonts();
	pdf(file=paste0(wd,f_pheno_discrimination, f_results,f_paralellePlot, name,"_",nameOutput,".pdf"),height=height,width=width)
	print(plot_noScaling.tendency)
	print(plot_noScaling_points.tendency)
	print(plot_Standard.tendency)
	print(plot_StandardCenter.tendency)
	print(plot_Univariately.tendency)
	print(plot_noScaling.ind)
	print(plot_Standard.ind)
	print(plot_StandardCenter.ind)
	print(plot_noScaling.tendency.CI)
	print(plot_Standard.tendency.CI)
	print(plot_StandardCenter.tendency.CI)
	print(plot_Univariately.tendency.CI)
	dev.off()
	

	assign(paste0("plot_noScaling.tendency_",nameOutput),plot_noScaling.tendency,.GlobalEnv);
	assign(paste0("plot_noScaling_points.tendency_",nameOutput),plot_noScaling_points.tendency,.GlobalEnv);
	assign(paste0("plot_noScaling.ind_",nameOutput),plot_noScaling.ind,.GlobalEnv);
	assign(paste0("plot_noScaling.tendency.CI_",nameOutput),plot_noScaling.tendency.CI,.GlobalEnv);

	assign(paste0("plot_Standard.tendency_",nameOutput),plot_Standard.tendency,.GlobalEnv);
	assign(paste0("plot_Standard.ind_",nameOutput),plot_Standard.ind,.GlobalEnv);
	assign(paste0("plot_Standard.tendency.CI_",nameOutput),plot_Standard.tendency.CI,.GlobalEnv);

	assign(paste0("plot_StandardCenter.tendency_",nameOutput),plot_StandardCenter.tendency,.GlobalEnv);
	assign(paste0("plot_StandardCenter.ind_",nameOutput),plot_StandardCenter.ind,.GlobalEnv);
	assign(paste0("plot_StandardCenter.tendency.CI_",nameOutput),plot_StandardCenter.tendency.CI,.GlobalEnv);

	assign(paste0("plot_Univariately.tendency_",nameOutput),plot_Univariately.tendency,.GlobalEnv);
	assign(paste0("plot_Univariately.tendency.CI_",nameOutput),plot_Univariately.tendency.CI,.GlobalEnv);
	
};


