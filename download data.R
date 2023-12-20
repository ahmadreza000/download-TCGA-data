library(TCGAbiolinks)

#GDCProjects = as.data.frame(getGDCprojects())
#TCGA_projects = GDCProjects["CPTAC" == substr(GDCProjects$project_id,1,5),]

query <- GDCquery(
  project = "TCGA-THCA", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  sample.type =c('Primary Tumor',"Solid Tissue Normal")
)
resultOfData0 = getResults(query)
#table(resultOfData$sample_type)
GDCdownload(query)
write.csv(resultOfData0,"geo-THCA.csv")

query <- GDCquery(
  project = "TCGA-THCA", 
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification",
  sample.type =c("Primary Tumor",'Solid Tissue Normal')
)
resultOfData1 = getResults(query)
#table(resultOfData$sample_type)
GDCdownload(query)
write.csv(resultOfData1,"prot-THCA.csv")


query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "DNA Methylation",
  data.type = "Masked Intensities",
  platform = "Illumina Human Methylation 27",
  sample.type =c("Primary Tumor",'Solid Tissue Normal')
)
resultOfData2 = getResults(query)
#table(resultOfData$sample_type)
GDCdownload(query)
write.csv(resultOfData2,"meth27-LUAD.csv")


query.isoform <- GDCquery(
  project = "TCGA-THCA", 
  experimental.strategy = "miRNA-Seq",
  data.category = "Transcriptome Profiling", 
  data.type = "Isoform Expression Quantification",
  sample.type =c("Primary Tumor",'Solid Tissue Normal')
)

resultOfData3 = getResults(query.isoform)
#table(resultOfData$sample_type)
GDCdownload(query.isoform)
write.csv(resultOfData3,"miRNA-THCA.csv")


query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "DNA Methylation",
  data.type = "Masked Intensities",
  platform = "Illumina Human Methylation 450",
  sample.type =c("Primary Tumor",'Solid Tissue Normal')
)

resultOfData4 = getResults(query)
GDCdownload(query)
write.csv(resultOfData4,"met450-THCA.csv")
#table(resultOfData$sample_type)





query <- GDCquery(
  project = "TCGA-LGG", 
  data.category = "Clinical", 
  data.format = "bcr xml",
)
GDCdownload(query)
clinical_followup <- GDCprepare_clinic(query,clinical.info=c("follow_up"))
clinical_patiant <- GDCprepare_clinic(query,clinical.info=c("patient"))
clinical_stage_event <- GDCprepare_clinic(query,clinical.info="stage_event")
clinical_drug<- GDCprepare_clinic(query,clinical.info="drug")
clinical_new_tumor_event<- GDCprepare_clinic(query,clinical.info="new_tumor_event")
clinical <- GDCquery_clinic(project = "TCGA-LGG", type = "clinical")

write.csv(clinical,"clinical-LGG.csv",row.names = F)
write.csv(clinical_followup,"clinical_followup-LGG.csv",row.names = F)
write.csv(clinical_patiant,"clinical_patiant-LGG.csv",row.names = F)
write.csv(clinical_drug,"clinical_drug-LGG.csv",row.names = F)