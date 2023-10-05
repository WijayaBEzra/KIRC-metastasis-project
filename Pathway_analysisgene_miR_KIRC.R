# Pathway analysis
#
# given the genes set significantlly affected by miRNA we'd like to check its impact in both N-P, and P-M based on gene expression analysis
# there are several tools : https://cran.r-project.org/web/packages/pathfindR/vignettes/intro_vignette.html
# for manual usege : https://cran.r-project.org/web/packages/pathfindR/vignettes/manual_execution.html
# OR library(pathview); Note:  Particullary, users are required to formally cite the original Pathview paper
# tools mapping its pathway : https://rpubs.com/barryus/class15 

#Start : Gene Expression analysis of selected gene [those who have significant correlated Target genes miRNA only]

###Required Installation###
#install.packages("pathfindR")
#BiocManager::install("KEGGREST")
#BiocManager::install("KEGGgraph")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("org.Hs.eg.db")

#BiocManager::install("pathview")
#BiocManager::install("gage")
#BiocManager::install("gageData")

library(pathfindR)
#error: Java version not detected. Please download and install Java from “https://www.java.com/en/”
#After installation of JAVE (latest) pleaser restart your R
library(pathview) #please cite his paper
library(gage)
library(gageData)

library("AnnotationDbi")
library("org.Hs.eg.db")

library(dplyr)
library(stringr)
library(tibble)

library(limma)
library(EnhancedVolcano)

options(scipen = 9)

### --------------------------------------- DATA EXTRACTION Gene Expression -------------------------------------------###
# NOTE : data value is log2(count+1)
# Only choose read_count criteria
myFile = read.delim(file.choose(), header = TRUE) #WARNING: check your tab, choose your file manually to Gene Expression
Gene_name=data.frame(myFile[,1])

### Divide Data into normal - Primary - Metastasis ###
#Seperate Normal-Primary Cancer
Normal_gene=myFile[,grep(".11", colnames(myFile))] #miRNA name is lost in this data
Primary_gene=myFile[,grep(".01", colnames(myFile))] #need to grep patient of metastasis name from clinical

#Get Clinical data metastasis
myFile2 = read.delim(file.choose(), header = TRUE) #warning: check your tab, choose your clinical file manually
new_row= data.frame(myFile2$sampleID,myFile2$pathologic_M) #only choose col SampleID & TNM (M) stage
cc=dplyr::filter(new_row, grepl('M1', myFile2$pathologic_M)) #Grab row only M1
cc_=str_replace_all(cc$myFile2.sampleID,"-",".") # edited TCGA name in from clinical data which use '-' in the patient's name
da=(paste(cc_, collapse="|")) #To make it as "OR"
Metastasis_gene=Primary_gene %>% dplyr::select(matches(da)) #calling metastasis data from Primary Data
row.names(Metastasis_gene)=Gene_name[,1]

#Update Primary M0 only
cp=dplyr::filter(new_row, grepl('M0', myFile2$pathologic_M)) #Grab patients row only M0
cp_=str_replace_all(cp$myFile2.sampleID,"-",".")
dap=(paste(cp_, collapse="|")) #To make it as "OR"
Primary_gene=Primary_gene %>% dplyr::select(matches(dap)) #take Patients name in clinical with M0 cathegory
row.names(Primary_gene)=Gene_name[,1]

### Limma analysis of DEGs ###
## DEG's Analysis P-N & M-P (data preparation) ##
design_n=data.frame(colnames(Normal_gene),code=0) ; colnames(design_n)=c("patient_ID","code")
design_p=data.frame(colnames(Primary_gene),code=1) ; colnames(design_p)=c("patient_ID","code")
design_m=data.frame(colnames(Metastasis_gene),code=2); colnames(design_m)=c("patient_ID","code")

# Design Matrix
design_pvn=rbind(design_p,design_n)
design_pvn=data.matrix(design_pvn, rownames.force = NA) #convert DF to numeric matrix
design_mvp=rbind(design_m, design_p)
design_mvp=data.matrix(design_mvp, rownames.force = NA) #convert DF to numeric matrix

# Preparing datasets combined data.frame
# Primary vs Normal
    df_primary_v_normal=data.frame(Gene_name,Primary_gene,Normal_gene)
    rownames(df_primary_v_normal) <- df_primary_v_normal[,1] #col1 as rownames
    df_primary_v_normal = df_primary_v_normal[,-1]
    df_primary_v_normal=(2^(df_primary_v_normal))-1
    df_primary_v_normal=data.matrix(df_primary_v_normal,rownames.force = NA) #Preparing value as numeric
    
    y <- voom(df_primary_v_normal, design_pvn, plot = T)
    dim(df_primary_v_normal)
    
    # DEG Analaysis
    fit <- lmFit(y, design_pvn)
    fit <- eBayes(fit)
    topTable(fit)
    
    top.table <- topTable(fit, sort.by = "F", n = Inf)
    
    # Visual simple volcanoplot
    volcanoplot(fit, coef=2, names = fit$F, 
                highlight = 5)
    
    toptable <- topTable(fit, n = Inf)
    toptable=add_column(toptable, row.names(toptable), .after = 0) #adding gene name in toptable
    
    #more filter on volcanoplot
    # link : https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
    # Warning message:
    # One or more p-values is 0. Converting to 10^-1 * current lowest non-zero p-value... 
    EnhancedVolcano(toptable,
                    lab = row.names(toptable),
                    x = 'code',
                    y = 'P.Value',
                    title = 'Primary vs Normal Gene Expression',
                    pCutoff = 10e-32,
                    FCcutoff = 0.5,
                    pointSize = 3.0, #dot font size
                    labSize = 3.0) #mirna name font size
    
    #fixed the table names
    toptable_=toptable
    colnames(toptable_)=c("Gene", "patient_ID",  "logFC",       "AveExpr",    "F" ,         "P.Value",    "adj.P.Val" )
    Result_pvn=data.frame(toptable_) #Table results

# Metastasis vs Primary 
    df_metastasis_v_primary=data.frame(Gene_name, Metastasis_gene, Primary_gene)
    rownames(df_metastasis_v_primary) <- df_metastasis_v_primary[,1] #col1 as rownames
    df_metastasis_v_primary = df_metastasis_v_primary[,-1]
    df_metastasis_v_primary=(2^(df_metastasis_v_primary))-1
    df_metastasis_v_primary=data.matrix(df_metastasis_v_primary,rownames.force = NA) #Preparing value as numeric
    
    x <- voom(df_metastasis_v_primary, design_mvp, plot = T)
    dim(df_metastasis_v_primary)
    
    # DEG Analaysis
    fitm <- lmFit(x, design_mvp)
    fitm <- eBayes(fitm)
    topTable(fitm)
    
    top.table_m <- topTable(fitm, sort.by = "F", n = Inf)
    
    # Visual simple volcanoplot
    volcanoplot(fitm, coef=2, names = fitm$F, 
                highlight = 5)
    
    toptable_m <- topTable(fitm, n = Inf)
    toptable_m=add_column(toptable_m, row.names(toptable_m), .after = 0) #adding gene name in toptable
    
    #more filter on volcanoplot
    # link : https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
    # Note some p-value are too small which equals to zero
    # HIDE zero status p-value & adj pvalue
    toptable_mf=toptable_m
    toptable_m_filtered <- filter(toptable_mf, toptable_mf$adj.P.Val==0)# takes only zero adj.pval
    toptable_m_filtered2 <- filter(toptable_mf, toptable_mf$adj.P.Val>0) #takes only pval>0
    
    EnhancedVolcano(toptable_m_filtered2,
                    lab = row.names(toptable_m_filtered2),
                    x = 'code',
                    y = 'P.Value',
                    title = 'Metastasis vs Primary Gene Expression',
                    pCutoff = 10e-32,
                    FCcutoff = 0.5,
                    pointSize = 3.0, #dot font size
                    labSize = 3.0) #mirna name font size
    
    #fixed the table names
    toptable_m_=toptable_m
    colnames(toptable_m_)=c("Gene", "patient_ID",  "logFC",       "AveExpr",    "F" ,         "P.Value",    "adj.P.Val" )
    Result_mvp=data.frame(toptable_m_) #Table results

#Quick SAVE results
library(writexl)
write_xlsx(Result_pvn,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Result_pvn_geneexpression.xlsx")
write_xlsx(Result_mvp,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Result_mvp_geneexpression.xlsx")    

### --------------------------------------- Get miRNA-GENE Sig Correlation -------------------------------------------###
# note : Every miRNA has unique correlation for each genes BUT in this case we combined all target genes in each miRNA to see what are those gene affecting in biological function pathway
# Load Previous Results 
library(openxlsx) #install these package if not exist in your R

#get primary mir
sheets <- openxlsx::getSheetNames("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Primary_Filtered_correlation_mir_tg.xlsx")
primary_mir <- lapply(sheets, openxlsx::read.xlsx, xlsxFile="D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Primary_Filtered_correlation_mir_tg.xlsx")
names(primary_mir) <- sheets
# grab specifics miRNA from data lists
mir19a_tg=data.frame(primary_mir$mir19a$row.names.cor_mir19a_tg); colnames(mir19a_tg)=c("target_gene") #get all col miRNA correlate target genes
mir101_tg=data.frame(primary_mir$mir101$row.names.cor_mir101_tg.); colnames(mir101_tg)=c("target_gene") #get tg mirna 101
mir3065_tg=data.frame(primary_mir$mir3065$row.names.cor_mir3065_tg.); colnames(mir3065_tg)=c("target_gene") #get tg mirna 3065
mir3130_tg=data.frame(primary_mir$mir3130$row.names.cor_mir3130_tg.); colnames(mir3130_tg)=c("target_gene") #get tg miRNA 3130

miRNA_tg_primary = rbind(mir19a_tg,mir101_tg,mir3065_tg,mir3130_tg) #combined all mirna target genes
# QC check if there are duplicate
miRNA_tg_primary[!duplicated(miRNA_tg_primary), ]


#get metas mir
sheets_m <- openxlsx::getSheetNames("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Metastasis_Filtered_correlation_mir_tg.xlsx")
metas_mir <- lapply(sheets_m, openxlsx::read.xlsx, xlsxFile="D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Metastasis_Filtered_correlation_mir_tg.xlsx")
names(metas_mir) <- sheets_m

mir10b_tg=data.frame(metas_mir$mir10b$row.names.cor_mir10b_tg.); colnames(mir10b_tg)=c("target_gene")
mir99b_tg=data.frame(metas_mir$mir99b$row.names.cor_mir99b_tg.); colnames(mir99b_tg)=c("target_gene")
mir21_tg=data.frame(metas_mir$mir21$row.names.cor_mir21_tg.); colnames(mir21_tg)=c("target_gene")
mir149_tg=data.frame(metas_mir$mir149$row.names.cor_mir149_tg.); colnames(mir149_tg)=c("target_gene")
mir410_tg=data.frame(metas_mir$mir410$row.names.cor_mir410_tg.); colnames(mir410_tg)=c("target_gene")
mir223_tg=data.frame(metas_mir$mir223$row.names.cor_mir223_tg.); colnames(mir223_tg)=c("target_gene")
mir4746_tg=data.frame(metas_mir$mir4746$row.names.cor_mir4746_tg.); colnames(mir4746_tg)=c("target_gene")
mir491_tg=data.frame(metas_mir$mir491$row.names.cor_mir491_tg.); colnames(mir491_tg)=c("target_gene")
mir200b_tg=data.frame(metas_mir$mir200b$row.names.cor_mir200b_tg.); colnames(mir200b_tg)=c("target_gene")
mir92b_tg=data.frame(metas_mir$mir92b$row.names.cor_mir92b_tg.); colnames(mir92b_tg)=c("target_gene")

miRNA_tg_metastasis = rbind(mir10b_tg,mir99b_tg,mir21_tg,mir149_tg,mir410_tg,mir223_tg,mir4746_tg,mir491_tg,mir200b_tg,mir92b_tg) #combined all mirna target genes

#Prepare for data filtration for DEG's
pa=t(miRNA_tg_primary)
pa=(paste(pa, collapse="|")) #To make it as "OR"
ma=t(miRNA_tg_metastasis)
ma=(paste(ma, collapse="|")) #To make it as "OR"

pvn=t(Result_pvn) #transpose so gene as col
mvp=t(Result_mvp) #transpose so gene as col
Primary_target_gene=data.frame(pvn) %>% dplyr::select(matches(pa)) # col - col filter 
Metastasis_target_gene=data.frame(mvp) %>% dplyr::select(matches(ma)) # col - col filter 

Primary_target_gene=t(Primary_target_gene)
Metastasis_target_gene=t(Metastasis_target_gene)

#QC after transpose it turn into char instead of data frame
Primary_target_gene=data.frame(Primary_target_gene)
Metastasis_target_gene=data.frame(Metastasis_target_gene)

### ---------------------------------- Pathway analysis of miRNA-TARGET GENE Sig Correlation -----------------------------------------###
## Active Subnetwork Search and Enrichment Analyses ##
Primary_target_gene_=data.frame(Primary_target_gene$Gene,Primary_target_gene$logFC,Primary_target_gene$adj.P.Val )
Metastasis_target_gene_=data.frame(Metastasis_target_gene$Gene, Metastasis_target_gene$logFC,Metastasis_target_gene$adj.P.Val)

# Input col needs to be [OLD GENE | GENE | CHANGE | P_VALUE]
colnames(Primary_target_gene_)=c("GENE","CHANGE","P_VALUE")
colnames(Metastasis_target_gene_)=c("GENE","CHANGE","P_VALUE")


#Try another code from (https://cran.r-project.org/web/packages/pathfindR/vignettes/intro_vignette.html)
#   output_df <- run_pathfindR(Primary_target_gene_, gene_sets = "KEGG")
#error Error in pathfindR::input_testing(input, p_val_threshold) : p values must all be numeric

Primary_target_gene_a <- Primary_target_gene_ %>% mutate_at(c('CHANGE', 'P_VALUE'), as.numeric) #set certain col as numeric
output_df <- run_pathfindR(Primary_target_gene_a, gene_sets = "KEGG") #NEED TO WAIT +/ 5 min
#Note: this code although takes time but it provides up to Picture of Pathway

Metastasis_target_gene_a <- Metastasis_target_gene_ %>% mutate_at(c('CHANGE', 'P_VALUE'), as.numeric) #set certain col as numeric
output_df_m <- run_pathfindR(Metastasis_target_gene_a, gene_sets = "KEGG") #NEED TO WAIT +/ 5 min
#Note: this code although takes time but it provides up to Picture of Pathway

visualize_terms(result_df = output_df, 
                hsa_KEGG = FALSE, # boolean to indicate whether human KEGG gene sets were used for enrichment analysis or not
                pin_name_path = "KEGG")
#note : Stored in R workspaces [C:\Users\Ezra Bernardus\Documents\term_visualizations] -all related pathway

visualize_terms(result_df = output_df_m, 
                hsa_KEGG = FALSE, # boolean to indicate whether human KEGG gene sets were used for enrichment analysis or not
                pin_name_path = "KEGG")

# Pathway comparison between Primary & Metastasis
# Vector of "Case" IDs for Primary
cases <- colnames(Metastasis_gene)

#Filtered needed genes BASED on METASTASIS TARGET GENE PATHWAY
MvP_exp_mat=t(df_metastasis_v_primary) #Transpose Matrix
Metas_tg_filter=t(Metastasis_target_gene[,1])
Metas_tg_filter=(paste(Metas_tg_filter, collapse="|")) #To make it as "OR"
MvP_exp_mat=data.frame(MvP_exp_mat) %>% dplyr::select(matches(Metas_tg_filter)) # col - col filter 
MvP_exp_mat=t(MvP_exp_mat) #Return to Col=patients name; row = as gene

## Calculate scores for representative terms 
## PLOT HEATMAP using term descriptions
M_clustered <- cluster_enriched_terms(output_df_m, plot_dend = FALSE, plot_clusters_graph = FALSE) # Cluster results
#Note : The maximum average silhouette width was 0.17 for k = 31 
score_matrix <- score_terms(enrichment_table = M_clustered[M_clustered$Status == "Representative", ], #need to check status
                            exp_mat = MvP_exp_mat, #ALL expression sample both Metastasis-Primary [as.number]
                            cases = cases,        #Patients name in  Category
                            use_description = TRUE, # default FALSE
                            label_samples = FALSE, # default = TRUE
                            case_title = "Metastasis",  # default = "Case"
                            control_title = "Primary", # default = "Control"
                            low = "#f7797d", # default = "green"
                            mid = "#fffde4", # default = "black"
                            high = "#1f4037") # default = "red"
# Note: value is in count; Cluster no color in Plots view because too many patients to show, if zoom we can see the diff each patients 

### For Verification only 17.04.2023 ###
#
#
# we check unique set of Metastasis & Primary
# note:  using https://bioinformatics.psb.ugent.be/webtools/Venn/ We find set and union of pathway enriched PvN with MvP
# But in this case we want to check SETs in Metastasis only; Input was from Venn web-aplication

see_a=M_clustered[M_clustered$Status =="Representative",]
Set_onlyMetastasis = c("hsa04810",	"hsa00280",	"hsa04261",	"hsa05030",	"hsa00533",	"hsa05202",
                       "hsa04216",	"hsa04350",	"hsa04062",	"hsa04114",	"hsa04620",	"hsa04923",
                       "hsa04925",	"hsa04918",	"hsa05145",	"hsa04217",	"hsa04971",	"hsa04137",
                       "hsa04514",	"hsa04330",	"hsa04928",	"hsa05133",	"hsa04927",	"hsa05169",
                       "hsa03030",	"hsa04612",	"hsa04933",	"hsa04215",	"hsa04611",	"hsa05207",
                       "hsa04659",	"hsa04360",	"hsa04931",	"hsa04510",	"hsa04211",	"hsa04020",
                       "hsa04152" ,	"hsa00270",	"hsa04130",	"hsa04066",	"hsa04390",	"hsa04115",
                       "hsa04962",	"hsa03022",	"hsa04970",	"hsa03410",	"hsa05142",		
                       "hsa04022",	"hsa04520",	"hsa04024",	"hsa05171",	"hsa00970",		
                       "hsa05166",	"hsa05135",	"hsa04924",	"hsa04730",				
                       "hsa04392")								

Filtered_M_clustered = subset(M_clustered,grepl(paste0(Set_onlyMetastasis, collapse = "|"), M_clustered$ID))
score_matrix_fil <- score_terms(enrichment_table = Filtered_M_clustered, #need to check status for Representative & pathway member
                            exp_mat = MvP_exp_mat, #ALL expression sample both Metastasis-Primary [as.number]
                            cases = cases,        #Patients name in  Category in here is metastasis
                            use_description = TRUE, # default FALSE
                            label_samples = FALSE, # default = TRUE
                            case_title = "Metastasis",  # default = "Case"
                            control_title = "Primary", # default = "Control"
                            low = "#f7797d", # default = "green"
                            mid = "#fffde4", # default = "black"
                            high = "#1f4037") # default = "red"
# Note : results not satifying [removing weak value in data frame] 
# Try build treshold 
#try log each datasest in matrix
MvP_exp_mat_log2=log2((MvP_exp_mat)+1) #option-1 log method
MvP_exp_mat_scale <- as.data.frame(scale(MvP_exp_mat)) #standar scaling where mean will be ZERO values x<mean+SD Or x>mean+SD 
    MvP_exp_mat_scale  = as.matrix(MvP_exp_mat_scale )

#randomized patienst for Primary & metastasis
statistics_M=MvP_exp_mat_scale[,1:79] #Metastasis patients only
statistics_M_rand=statistics_M[,sample(ncol(statistics_M), size=50)] #randomly choose col from Metastasis patients
cases_rand=colnames(statistics_M_rand) #adjusting metastasis case accoridng to random pick

statistics_P=MvP_exp_mat_scale[,80:502] #Primary only
statistics_P_rand=statistics_P[,sample(ncol(statistics_P), size=50)] #randomly choose col from primary patients

    statistics_rand_model=cbind(statistics_M_rand,statistics_P_rand)
    
    
score_matrix_fil2 <- score_terms(enrichment_table = Filtered_M_clustered[Filtered_M_clustered$Status == "Representative", ], #need to check status for Representative & pathway member
                                exp_mat = statistics_rand_model, #ALL expression sample both Metastasis-Primary [as.number]
                                cases = cases,        #Patients name in  Category in here is metastasis
                                use_description = TRUE, # default FALSE
                                label_samples = FALSE, # default = TRUE
                                case_title = "Metastasis",  # default = "Case"
                                control_title = "Primary", # default = "Control"
                                low = "#f7797d", # default = "green"
                                mid = "#fffde4", # default = "black"
                                high = "#1f4037") # default = "red"

score_matrix_fil2 <- score_terms(enrichment_table = Filtered_M_clustered, #ALL Representative & pathway member
                                 exp_mat = statistics_rand_model, #ALL expression sample both Metastasis-Primary [as.number]
                                 cases = cases_rand,        #Patients name in  Category in here is metastasis
                                 use_description = TRUE, # default FALSE
                                 label_samples = FALSE, # default = TRUE
                                 case_title = "Metastasis",  # default = "Case"
                                 control_title = "Primary", # default = "Control"
                                 low = "#f7797d", # default = "green"
                                 mid = "#fffde4", # default = "black"
                                 high = "#1f4037") # default = "red"

#Focus : Based on Journal Review Core Metastasis Pathway
Set_JReview_Metastasis = c("hsa04510", "hsa04514", "hsa04810", "hsa04064", "hsa04310", 
                          "hsa04064", "hsa04012", "hsa04015", "hsa04020", "hsa04620", 
                          "hsa04621", "hsa04622", "hsa04010", "hsa04151","hsa04330", "hsa04310")		#update 16.05.2023 add PI3K						
Fil_JReview_Metastasis= subset(M_clustered,grepl(paste0(Set_JReview_Metastasis, collapse = "|"), M_clustered$ID))

# Random Log2 Normalization
statistics_M_log=MvP_exp_mat_log2[,1:79] #Metastasis patients only
statistics_M_rand_log=statistics_M_log[,sample(ncol(statistics_M_log), size=50)] #randomly choose col from Metastasis patients
cases_rand_log=colnames(statistics_M_rand_log) #adjusting metastasis case accoridng to random pick
statistics_P_log=MvP_exp_mat_log2[,80:502] #Primary only
statistics_P_rand_log=statistics_P_log[,sample(ncol(statistics_P_log), size=50)] #randomly choose col from primary patients
statistics_rand_model_log=cbind(statistics_M_rand_log,statistics_P_rand_log)


# using normalization MEAN
score_matrix_fil3 <- score_terms(enrichment_table = Fil_JReview_Metastasis, #ALL Representative & pathway member
                                 exp_mat = statistics_rand_model, #ALL expression sample both Metastasis-Primary [as.number]
                                 cases = cases_rand,        #Patients name in  Category in here is metastasis
                                 use_description = TRUE, # default FALSE
                                 label_samples = FALSE, # default = TRUE
                                 case_title = "Metastasis",  # default = "Case"
                                 control_title = "Primary", # default = "Control"
                                 low = "#f7797d", # default = "green"
                                 mid = "#fffde4", # default = "black"
                                 high = "#1f4037") # default = "red"

# using normalization LOG2
score_matrix_fil3 <- score_terms(enrichment_table = Fil_JReview_Metastasis, #ALL Representative & pathway member
                                 exp_mat = statistics_rand_model_log, #ALL expression sample both Metastasis-Primary [as.number]
                                 cases = cases_rand_log,        #Patients name in  Category in here is metastasis
                                 use_description = TRUE, # default FALSE
                                 label_samples = FALSE, # default = TRUE
                                 case_title = "Metastasis",  # default = "Case"
                                 control_title = "Primary", # default = "Control"
                                 low = "#f7797d", # default = "green"
                                 mid = "#fffde4", # default = "black"
                                 high = "#1f4037") # default = "red"

### ---------------------------------- Pathway visualization KEGG -----------------------------------------###
# Required Fold changes data [note: we only use metastasis because we'd like to compare primary vs metastasis]
foldchanges=Metastasis_target_gene$logFC
foldchanges=data.frame(foldchanges)

# ** GENE NAME CONVERSION CODE
Entrenz_genename = mapIds(org.Hs.eg.db,
                            keys=Metastasis_target_gene$Gene , 
                            column="ENTREZID",
                            keytype="SYMBOL",
                            multiVals="first")
Entrenz_genename_=data.frame(Entrenz_genename)
#NOTE : Some gene names is missing NA (which maybe duplicates)

foldchangesx = data.frame(Entrenz_genename_,foldchanges) #combined gene name (entrez) and their logFC
foldchangesx_ <- foldchangesx[rowSums(is.na(foldchangesx)) == 0, ] # Remove NA row in logFC metastasis data
row.names(foldchangesx_)=foldchangesx_$Entrenz_genename # switch to rownames
foldchangesx_[1] = NULL #use data from RA expreiment
foldchangesx_ <- foldchangesx_ %>% mutate_at(c('foldchanges'), as.numeric) #set logfold col as numeric
head(foldchangesx_)

data(kegg.sets.hs) #load kegg datasets

# Which will give us separate lists for pathways that are upregulated versus pathways that are down-regulated.
keggres = gage(foldchangesx_, gsets=kegg.sets.hs, same.dir=TRUE) #keggres accept entrez name only
attributes(keggres) # [1] "greater" "less"    "stats"
head(keggres$greater)

# show enrichment score
Up_pathway_kegg=data.frame(row.names(keggres$greater), keggres$greater)
Down_pathway_kegg=data.frame(row.names(keggres$less),keggres$less)

# Visualized alteration in Pathway
pathview(gene.data=foldchangesx_, pathway.id="hsa05211") #hsa05211 is about "Renal cell carcinoma"
#stored in Documents AS JPG

### ----------------------------------------SAVE all works--------------------- ###
library(writexl)
# Results -library(pathfindR): focus on overall pathway in logFC
write_xlsx(output_df,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\pathwayresult_DEG_TG_MIR_PrimaryVNormal.xlsx")
write_xlsx(output_df_m,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\pathwayresult_DEG_TG_MIR_MetastasisVPrimary.xlsx")

# Results - library(pathview) : focus on visualisation in KEGG
write_xlsx(Up_pathway_kegg,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\pathwayresult_UP_MetasKEGG.xlsx" )
write_xlsx(Down_pathway_kegg,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\pathwayresult_Down_MetasKEGG.xlsx" )

# additional save : 
# 1. heatmap in directory : D:\Dissertation 4.08.2022\Publication - 1 Bioinformatics miRNA
# 2. pathway with up & down gene color : Documents [R-directory]
# 3. Pathway graph using PathfindR stored in folder "term_visualizations"
# 4. similar result as table output_df stored in folder document : "pathfindR_Results(2)"


































### ------------[DOESNT WORK]---------------------- Pathway analysis of miRNA-TARGET GENE Sig Correlation -----------------------------------------###
## set input data according to pathfindR 
Primary_target_gene_=data.frame(Primary_target_gene$Gene,Primary_target_gene$logFC,Primary_target_gene$adj.P.Val )
Metastasis_target_gene_=data.frame(Metastasis_target_gene$Gene, Metastasis_target_gene$logFC,Metastasis_target_gene$adj.P.Val)
colnames(Primary_target_gene_)=c("Gene.symbol","logFC ","adj.Pval")
colnames(Metastasis_target_gene_)=c("Gene.symbol","logFC ","adj.Pval")

## Filter only significant genes [wont be needed since QC: all is significant]
## filter gene used in biogrid network ; more about biogrid : https://thebiogrid.org/

Path_Primary_target_gene <- input_processing(input = Primary_target_gene_, # the input: in this case, differential expression results
                                 p_val_threshold = 0.05, # p value threshold to filter significant genes
                                 pin_name_path  = "STRING", # the name of the PIN to use for active subnetwork search
                                 #other Option : “Biogrid”, “STRING”, “GeneMania”, “IntAct”, “KEGG”, “mmu_STRING”
                                 convert2alias = TRUE) # boolean indicating whether or not to convert missing symbols to alias symbols in the PIN
#ERROR : pathfindR cannot handle p values < 1e-13. These were changed to 1e-13, None of the genes were in the PIN [solution: set data frame according required data frame, colomn title needs to be exactly the same]

#Using "BioCarta" as our gene sets for enrichment Analysis
# NOte in package : The available gene sets in pathfindR are “KEGG”, “Reactome”, “BioCarta”, “GO-All”, “GO-BP”, “GO-CC” and “GO-MF”
# more about Biocarta :https://maayanlab.cloud/Harmonizome/dataset/Biocarta+Pathways

#-----------------CHECK DATABASED---------------------#
biocarta_list <- fetch_gene_set(gene_sets = "BioCarta",
                                min_gset_size = 10,
                                max_gset_size = 300)
biocarta_gsets <- biocarta_list[[1]]
biocarta_descriptions <- biocarta_list[[2]]

# We Will try to use KEGG
KEGG_list <- fetch_gene_set(gene_sets = "KEGG",
                                min_gset_size = 10,
                                max_gset_size = 1000) #total 324 pathway listed in KEGG
KEGG_gsets <- KEGG_list[[1]]
KEGG_descriptions <- KEGG_list[[2]]

# We Will try to use Reactome
Reactome_list <- fetch_gene_set(gene_sets = "Reactome",
                            min_gset_size = 10,
                            max_gset_size = 1000) #total 324 pathway listed in Reactome
Reactome_gsets <- Reactome_list[[1]]
Reactome_descriptions <- Reactome_list[[2]]

# We Will try to use GO-ALL
GO_All_list <- fetch_gene_set(gene_sets = "GO-All",
                                min_gset_size = 10,
                                max_gset_size = 1000) #total 324 pathway listed in Reactome
GO_All_gsets <- GO_All_list[[1]]
GO_All_descriptions <- GO_All_list[[2]]

#------------------------------------------------------#

n_iter <- 2 ## number of iterations #adjust for trying
combined_res <- NULL ## to store the result of each iteration


## Active Subnetwork Search and Enrichment Analyses ##
# Input col needs to be [OLD GENE | GENE | CHANGE | P_VALUE]
colnames(Primary_target_gene_)=c("GENE","CHANGE","P_VALUE")
colnames(Metastasis_target_gene_)=c("GENE","CHANGE","P_VALUE")

# Active to find gene sets according to KEGG|Biocarta
for (i in 1:n_iter) {
  ###### Active Subnetwork Search
  snws_file <- paste0("active_snws_", i) # Name of output file
  active_snws <- active_snw_search(input_for_search = Path_Primary_target_gene, 
                                   pin_name_path = "STRING", 
                                   snws_file = snws_file,
                                   score_quan_thr = 0.8, # you may tweak these arguments for optimal filtering of subnetworks
                                   sig_gene_thr = 0.02, # you may tweak these arguments for optimal filtering of subnetworks
                                   search_method = "GR") # method of algorithm search : greedy search (GR)
  # All search algorithms: Options are the greedy search (GR), simulated annealing (SA), or genetic algorithm (GA)
  ###### Enrichment Analyses
  current_res <- enrichment_analyses(snws = active_snws,
                                     sig_genes_vec = Path_Primary_target_gene$GENE,
                                     pin_name_path = "STRING", 
                                     genes_by_term = GO_All_gsets,
                                     term_descriptions = GO_All_descriptions,
                                     adj_method = "bonferroni",
                                     enrichment_threshold = 0.05,
                                     list_active_snw_genes = TRUE) # listing the non-input active snw genes in output
  
  ###### Combine results via `rbind`
  combined_res <- rbind(combined_res, current_res)
}

# Summarize Combined Enrichment Results
summarized_df <- summarize_enrichment_results(combined_res, 
                                              list_active_snw_genes = TRUE)
#note: not find the active subnetworks???; TRY : filter according to path_primary_target_gene

#experiment 22.03.2023 :
# assuming pin_name path is : Biogrid
# -KEGG : found 0 active subnetworks
# -biocarta : found 0 active subnetworks???? why?
# -Reactome : found 0
# -GO-ALL : found 0

# assuming pin_name path is : STRING
# -KEGG : found 0 active subnetworks
# -biocarta : found 0 active subnetworks???? why?
# -Reactome : found 0
# -GO-ALL : found 0

# QC check using msigDB
#library(writexl)
#write_xlsx(Primary_target_gene_,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\test.xlsx")
# Manual checkcan find but not common target genes
