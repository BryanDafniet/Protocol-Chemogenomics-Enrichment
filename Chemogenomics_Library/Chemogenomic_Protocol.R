#library(ggplot2)
library(dplyr)
#library(viridis)
#library(hrbrthemes)
#library(RColorBrewer)
#library(ggdendro)
library(factoextra)
library(NbClust)
#library(cowplot)
library(stringr)

#Scaffold selection

# Select scaffold with most target and, if equal, the lower distance based on the distance matrix and the cluster it is in.
Most_target <- function(indiv) {
  SommeMax=0
  TargetMax=0
  IndivTarget=""
  Tmp_order=as.data.frame(apply(df.dist_uni[indiv,indiv],1,sum))
  #  print(df.dist_uni[indiv,indiv])
  colnames(Tmp_order)[1]="Score"
  Tmp_order[,1]=Tmp_order[order(Tmp_order[,1]),]
  All_best_score=rownames(Tmp_order)[Tmp_order[,1]==Tmp_order[1,1]]
  print(Tmp_order)
  print(All_best_score)
  for (individu in All_best_score) {
    tmp_comp=as.data.frame(ScaTarget[which(ScaTarget[,"Scaffold"]==individu),c(1,2,3)])
    print(tmp_comp)
    sum_tmp=sum(tmp_comp[,3])
    tmp_tar1=dim(tmp_comp)[1]
    if (tmp_tar1>TargetMax) {
      IndivTarget=individu
      TargetMax=tmp_tar1
      SommeMax=sum_tmp
    } else if (tmp_tar1==TargetMax) {
      if (sum_tmp>SommeMax) {
        SommeMax=sum_tmp
        IndivTarget=individu
      }
    } else {
      #next
    }
    
  }
  print(IndivTarget)
  return(IndivTarget)
}

#Neo4j command to extract TargetSca2
# MATCH (g:ProteinClass)<-[:memberOf*1..]-(d:ProteinClass)
# WHERE NOT (g)-[:memberOf]->(:ProteinClass)
# with d
# MATCH (d)<-[:memberOf*0..]-(p:ProteinClass)<-[v:memberOf]-(u:UniprotInter)<-[i:memberOf]-(t:Target)<-[a:actifConf9]-(m:Molecule)<-[:substructureOf*2]-(sca:Scaffolds)
# WITH t.Name as TargetName, sca.ScaID as Scaffold, count(DISTINCT m.InchiKey) as Molecule
# RETURN DISTINCT TargetName, Scaffold,Molecule
#save as TargetSca2

ScaTarget=read.csv("TargetSca2.csv", header=T)
ScaTarget=ScaTarget[order(ScaTarget$Scaffold),]

sca.tmp=data.frame(matrix(ncol = 1, nrow = 32868, unique(ScaTarget$Scaffold)))
mat.target = data.frame(matrix(ncol = 1403, nrow = 32868, 0))
colnames(mat.target)=unique(ScaTarget$TargetName)
colnames(sca.tmp)[1]="Scaffold"
mat.target=cbind(sca.tmp,mat.target)
count = 1

#Retire les scaffold avec cible > 6 et change les 0 de la matrice en 1 quand le scaffold interagit avec une cible
for (sca in as.factor(unique(ScaTarget$Scaffold))) {
  #print(sca)
  tmp=ScaTarget[which(ScaTarget[,"Scaffold"]==sca),]
  tmp=tmp[1]
  #print(tmp)
  if (dim(tmp)[1] > 6) {
    mat.target=mat.target[-which(mat.target[,"Scaffold"]==sca),]
    next
    
  }
  for (target in tmp[,1]){
    mat.target[which(mat.target[,"Scaffold"]==sca),target]=replace(mat.target[which(mat.target[,"Scaffold"]==sca),target],values = 1)
    
  }
  count=count+1
  print(count)
  
}


remove.var = which(apply(mat.target,2,sd) == 0)
mat.target2 = mat.target[,-remove.var]
filtre=as.data.frame(colnames(mat.target2))
filtre=as.data.frame(filtre[-1,])

#Hclust for scaffolds/distance matrix
mat_distRetry=dist(x = mat.target2[,2:1222],method="binary")

mat_hclust=hclust(mat_distRetry, method="complete")

data_5k=cutree(mat_hclust, k = 5000)
names(data_5k) = as.character(mat.target2[which(rownames(mat.target2)==as.integer(names(data_5k))),1])
# data_5k[1:50]
table(is.na(names(data_5k))) ## If true error
mat_centr_hclust5 = c("k","Scaffold")
for(k in 1:5000){ # For each cluster we take the one the most in the middle
  indiv = names(data_5k[which(data_5k == k)])
  # If only one scaf per group we take the scaf
  # if 2 scaf we take the one with the highest number of target
  # if more the more central
  if(length(indiv)==0){
    # print(paste("next",k))
    next
  } else if(length(indiv)==1){
    # print(paste("length1",k))
    centroid5 = indiv
  } else {

    centroid5 = indiv[which(indiv==Most_target(indiv))]
  }
  mat_centr_hclust5 = rbind(mat_centr_hclust5,c(k,centroid5))
}
colnames(mat_centr_hclust5) = mat_centr_hclust5[1,]
mat_centr_hclust5 = mat_centr_hclust5[-1,]
mat_centr_hclust_fac5=as.data.frame(mat_centr_hclust_fac5)
mat_centr_hclust_fac5[,2]=as.factor(mat_centr_hclust_fac5[,2])

Sum_all5k=semi_join(ScaTarget, mat_centr_hclust_fac5,"Scaffold")
#write.csv(Sum_all5k, "5000_Comp.csv")

#Neo4J Command to extract information on molecules and protein class from the 5k scaffold
# MATCH (g:ProteinClass)<-[:memberOf*1..]-(d:ProteinClass)
# WHERE NOT (g)-[:memberOf]->(:ProteinClass)
# with d
# LOAD CSV WITH HEADERS FROM 'file:///5000_Comp.csv' as row FIELDTERMINATOR ','
# MATCH (d)<-[:memberOf*0..]-(p:ProteinClass)<-[v:memberOf]-(u:UniprotInter)<-[i:memberOf]-(t:Target)<-[a:actifConf9]-(m:Molecule)<-[:substructureOf*2]-(sca:Scaffolds)
# WHERE sca.ScaID = row.Scaffold AND t.Name = row.TargetName
# WITH p.Name as ProtClass,u.Name as UniprotInter, t.Name as Target, sca.ScaID as Scaffold, m.InchiKey as Molecule
# RETURN ProtClass,UniprotInter,Target, Scaffold, Molecule

#Save the data to AllData_5000Redone.csv

Scaff5k_Redone=read.csv("AllData_5000Redone.csv", header=T)

#Creation of the binary matrix molecule/target
count = 0
mol.tmp5=data.frame(matrix(ncol = 1, nrow = 41620, unique(Scaff5k_Redone$Molecule)))
mat.mol5 = data.frame(matrix(ncol = 850, nrow = 41620, 0))
colnames(mat.mol5)=unique(Scaff5k_Redone$UniprotInter)
colnames(mol.tmp5)[1]="Molecule"
mat.mol5=cbind(mol.tmp5,mat.mol5)
table(is.na(mat.mol5))
for (mol in as.factor(unique(Scaff5k_Redone$Molecule))) {
  print(mol)
  tmp=Scaff5k_Redone[which(Scaff5k_Redone[,"Molecule"]==mol),]
  tmp=unique(tmp[2])
  for (target in tmp[,1]){
    mat.mol5[which(mat.mol5[,"Molecule"]==mol),target]=replace(mat.mol5[which(mat.mol5[,"Molecule"]==mol),target],values = 1)
  }
  count=count+1
  print(count)
}

mat_dist_Mol5=dist(x = mat.mol5[2:851],method="binary")
table(is.na(mat_dist_Mol5))

df.dist_mol5=as.matrix(mat_dist_Mol5, labels=TRUE)
colnames(df.dist_mol5) <- rownames(df.dist_mol5) <- mat.mol5[['Molecule']]

#write.csv(mat.mol5,"mat_mol41k_Redone.csv")
#write.csv(df.dist_mol5,"Dist_Mat_Mol41k_Redone.csv")

#Pareto multiobjective optimization
###### Analysis and pareto #####

#Collect pareto file after the pareto optimization

Pareto=read.csv("Res_5000_1000_600.csv", header=T, sep="\t")
#ParetoLog=read.csv("Res_5000_1000_600log.csv", sep="\t",header=T)

#The pareto has a genetic algorithm base so results may differ. It should be similar but the subset selected and the objectives maximizations
# might differ slightly for the molecule selection.

Pareto1=Pareto[which(Pareto[,"ParetoFront"]==1),]
ParetoSub43=Pareto1[which(Pareto1[,"SubsetID"]==43),]


#Exploring proteins with no main family link
# MATCH (p:ProteinClass)<-[v:memberOf]-(u:UniprotInter)<-[i:memberOf]-(t:Target)<-[a:actifConf9]-(m:Molecule)<-[:substructureOf*2]-(sca:Scaffolds)
# WHERE NOT (p)-[:memberOf]->(:ProteinClass) AND p.Name <> "Membrane receptor" AND p.Name <> "Enzyme" AND p.Name <> "Ion channel" AND p.Name <> "Auxiliary transport protein" AND p.Name <> "Transporter" AND p.Name <> "Epigenetic regulator" AND p.Name <> "Transcription factor"
# WITH p.Name as ProtClass,u.Name as UniprotInter, t.Name as Target, sca.ScaID as Scaffold, m.InchiKey as Molecule
# RETURN DISTINCT ProtClass,UniprotInter, Molecule, Scaffold

#Save the results as OtherProtClass_WithMol.csv

Other_Prot_Comp=read.csv("OtherProtClass_WithMol.csv", )
Other_Prot_Comp=as.data.frame(unique(Other_Prot_Comp[,c(1,2)]))

# Check scaffolds and Molecules that are linked already to other unexplored proteins in our 5k dataset
Comp_Mol_Fam=unique(semi_join(Other_Prot_Comp,ParetoSub43,"Molecule"))
Comp_Mol_Fam=as.data.frame(unique(Comp_Mol_Fam[,2]))

Comp_Scaf_Fam=semi_join(Other_Prot_Comp,ParetoSub43,"Scaffold")
Comp_Scaf_Fam=as.data.frame(unique(Comp_Scaf_Fam[,2]))

colnames(Comp_Mol_Fam)[1]="Target"
colnames(Comp_Scaf_Fam)[1]="Target"

#Target on which a scaffold in the pareto selection was linked to, but the molecule selected related to that scaffold has no relationship with it
#We did not consider the question further as we considered a scaffold enough here to represent the hit on a target
Comp_ScaMol_Fam=anti_join(Comp_Scaf_Fam, Comp_Mol_Fam,"Target")
Comp_ScaMol_Fam_Target=as.data.frame(unique(Comp_ScaMol_Fam[,2]))
Comp_ScaMol_Fam=Comp_ScaMol_Fam[,-1]
Comp_ScaMol_Fam=unique(Comp_ScaMol_Fam)
#write.csv(Comp_ScaMol_Fam,"Scaffold_Switch.csv")

# unique(ProtToMap[,2])

## Check which target are not targeted by at least a scaffold in our pareto
ProtToMap=anti_join(Other_Prot_Comp,Comp_Scaf_Fam, "Target")
#PropToMap_Targetleft=unique(as.data.frame(ProtToMap[,2]))

#write.csv(ProtToMap,"ProtToMap.csv")


# LOAD CSV WITH HEADERS FROM 'file:///ProtToMap.csv' as row FIELDTERMINATOR ','
# MATCH (p:ProteinClass)<-[v:memberOf]-(u:UniprotInter)<-[i:memberOf]-(t:Target)<-[a:actifConf9]-(m:Molecule)<-[:substructureOf*2]-(sca:Scaffolds)
# WHERE NOT (p)-[:memberOf]->(:ProteinClass) AND p.Name <> "Membrane receptor" AND p.Name <> "Enzyme" AND p.Name <> "Ion channel" AND p.Name <> "Auxiliary transport protein" AND p.Name <> "Transporter" AND p.Name <> "Epigenetic regulator" AND p.Name <> "Transcription factor" AND u.Name = row.Target
# WITH p.Name as ProtClass,u.Name as Target, sca.ScaID as ScaID, m.InchiKey as Molecule
# RETURN DISTINCT Target,Molecule, ScaID
#save result as SummaryToMap.csv


AllInfo_ToMap=read.csv("SummaryToMap.csv", header=T)
AllInfo_ToMap=unique(AllInfo_ToMap[,c(2,3,4)])


#Using the SummaryToMap.csv on neo4j we will see if overhall our scaffold targets more than 6 proteins indicating their overspecificity
# LOAD CSV WITH HEADERS FROM 'file:///SummaryToMap.csv' as row FIELDTERMINATOR ','
# MATCH (p:ProteinClass)<-[v:memberOf]-(u:UniprotInter)<-[i:memberOf]-(t:Target)<-[a:actifConf9]-(m:Molecule)<-[:substructureOf*2]-(sca:Scaffolds)
# WHERE sca.ScaID= row.ScaID
# WITH u.Name as Target, sca.ScaID as ScaID, m.InchiKey as Molecule
# RETURN DISTINCT ScaID,Target, Count(DISTINCT Molecule)
#save as TotalTarget_Scaffold.csv

AllInfo_ToMap_Count=read.csv("TotalTarget_Scaffold.csv",header=T)
unique(AllInfo_ToMap_Count[,'ScaID'])
unique(AllInfo_ToMap_Count[,'Target'])
unique(AllInfo_ToMap[,'Molecule'])

sca.tmp_ToMap=data.frame(matrix(ncol = 1, nrow = 412, unique(AllInfo_ToMap_Count$ScaID)))
mat.target_ToMap = data.frame(matrix(ncol = 843, nrow = 412, 0))
colnames(mat.target_ToMap)=unique(AllInfo_ToMap_Count$Target)
colnames(sca.tmp_ToMap)[1]="Scaffold"
mat.target_ToMap=cbind(sca.tmp_ToMap,mat.target_ToMap)

#Retire les scaffold avec cible > 6 among the potential candidates to target the missing proteins
for (sca in as.factor(unique(AllInfo_ToMap_Count$ScaID))) {
  #print(sca)
  tmp=AllInfo_ToMap_Count[which(AllInfo_ToMap_Count[,"ScaID"]==sca),]
  tmp=tmp[1]
  #print(tmp)
  if (dim(tmp)[1] > 6) {
    mat.target_ToMap=mat.target_ToMap[-which(mat.target_ToMap[,"Scaffold"]==sca),]
    next
    
  }
}
mat.target_ToMap=as.data.frame(mat.target_ToMap[,1])
colnames(mat.target_ToMap)[1]="ScaID"
Check_Toadd=semi_join(AllInfo_ToMap,mat.target_ToMap,"ScaID" )
Check_NotTaken=anti_join(AllInfo_ToMap,Check_Toadd,"Target" )
#write.csv(Check_Toadd, "Activity_Toadd")

#The Scaffold_Switch.csv is from the Comp_ScaMol_Fam.
# LOAD CSV WITH HEADERS FROM 'file:///Scaffold_Switch.csv' as row FIELDTERMINATOR ','
# MATCH (sca:Scaffolds)-[:substructureOf*2]->(m:Molecule)-[:structureOf]->(c:CompoundName)-[]->(r:Result)-[:hasParameter]->(ap:AssayParameter)-[:parameterOf]->(a:Assay)-[:isOnTarget]->(t:Target)-[i:memberOf]->(u:UniprotInter)
# WHERE m.InchiKey = row.Molecule AND u.Name=row.Target AND a.Assay_type="B" AND a.ConfidenceScore="9" AND sca.ScaID=row.Scaffold
# WITH  u.Name as Target, m.InchiKey as Molecule, sca.ScaID as ScaID, ap.Parameter as Parameter, ap.Unit as Unit, r.Value as Value
# RETURN  DISTINCT Molecule,  Target, Value,  Unit,  Parameter, ScaID
#save as Switch_To_Select.csv

#Loop to select target with only one molecule and ,if multiple, select the one molecule with the highest activity
#The Activity_Mapped.csv is the same NEO4J command obtained from a the file Activity_Toadd.csv.
Activity_NM=read.csv("Activity_Mapped.csv", header = T) #The non mapped file
Test=anti_join(Activity_NM,Activity_Test,"Molecule")
unique(Activity_Test$Target)
unique(Activity_NM$Target)
Activity_NM=Activity_NM[which(Activity_NM[,"Parameter"] != "Inhibition"),]
Activity_IC50=Activity_NM[which(Activity_NM[,"Parameter"] == "IC50"),]
Activity_EC50=Activity_NM[which(Activity_NM[,"Parameter"] == "EC50"),]
Activity_Ki=Activity_NM[which(Activity_NM[,"Parameter"] == "Ki"),]


unique(Activity_NM[,2])
Select_OneMol = tibble()
for (target in as.factor(unique(Check_Toadd$Target))) {
  tmp=Check_Toadd[which(Check_Toadd[,"Target"]==target),]
  print(tmp)
  if (dim(tmp)[1]==1) {
    tmp=tmp[,c("Molecule", "Target", "ScaID")]
    colnames(tmp)[3]="ScaID"
    Select_OneMol=rbind(Select_OneMol, tmp[,c(2,1,3)])
  } else {
    To_split=tibble()
    tmp_sca=unique(as.data.frame(tmp[,3]))
    colnames(tmp_sca)[1]="ScaID"
    for(sca in tmp_sca$ScaID) {
      To_split=rbind(To_split,Activity_NM[which(Activity_NM[,"ScaID"]==sca & Activity_NM[,"Target"]==target),])
    }
    print(To_split)
    IC50=To_split[which(To_split[,"Parameter"]=="IC50"),]
    EC50=To_split[which(To_split[,"Parameter"]=="EC50"),]
    Ki=To_split[which(To_split[,"Parameter"]=="Ki"),]
    IC50=IC50[order(IC50[,"Value"]),]
    EC50=EC50[order(EC50[,"Value"]),]
    Ki=Ki[order(Ki[,"Value"]),]
    if (dim(IC50)[1]!=0) {
      print("IC50NOTNULL")
      Select_OneMol=rbind(Select_OneMol, IC50[1,c(2,1,6)])
    }
    if (dim(EC50)[1]!=0) {
      print("EC50NOTNULL")
      Select_OneMol=rbind(Select_OneMol, EC50[1,c(2,1,6)])
    }
    if (dim(Ki)[1]!=0) {
      print("KiNOTNULL")
      Select_OneMol=rbind(Select_OneMol, Ki[1,c(2,1,6)])
    }
  }
}


unique(Activity_NM[,2])

#From the targets that have a scaffold from our pareto subset linked to them
Select_OneMol=unique(Select_OneMol)

Activity=read.csv("Switch_To_Select.csv", header = T) # The switch file
Activity=Activity[which(Activity[,"Parameter"] != "Inhibition"),]
Activity_IC50=Activity[which(Activity[,"Parameter"] == "IC50"),]
Activity_EC50=Activity[which(Activity[,"Parameter"] == "EC50"),]
Activity_Ki=Activity[which(Activity[,"Parameter"] == "Ki"),]

unique(Activity[,2])
Select_OneMol_Switch = tibble()
for (target in as.factor(unique(Comp_ScaMol_Fam$Target))) {
  tmp=Comp_ScaMol_Fam[which(Comp_ScaMol_Fam[,"Target"]==target),]
  print(tmp)
  if (dim(tmp)[1]==1) {
    tmp=tmp[,c("Molecule", "Target", "Scaffold")]
    colnames(tmp)[3]="ScaID"
    Select_OneMol_Switch=rbind(Select_OneMol_Switch, tmp[,c(2,1,3)])
  } else {
    To_split=tibble()
    tmp_sca=unique(as.data.frame(tmp[,3]))
    colnames(tmp_sca)[1]="Scaffold"
    for(sca in tmp_sca$Scaffold) {
      To_split=rbind(To_split,Activity[which(Activity[,"ScaID"]==sca & Activity[,"Target"]==target),])
    }
    print(To_split)
    IC50=To_split[which(To_split[,"Parameter"]=="IC50"),]
    EC50=To_split[which(To_split[,"Parameter"]=="EC50"),]
    Ki=To_split[which(To_split[,"Parameter"]=="Ki"),]
    IC50=IC50[order(IC50[,"Value"]),]
    EC50=EC50[order(EC50[,"Value"]),]
    Ki=Ki[order(Ki[,"Value"]),]
    if (dim(IC50)[1]!=0) {
      print("IC50NOTNULL")
      Select_OneMol_Switch=rbind(Select_OneMol_Switch, IC50[1,c(2,1,6)])
    }
    if (dim(EC50)[1]!=0) {
      print("EC50NOTNULL")
      Select_OneMol_Switch=rbind(Select_OneMol_Switch, EC50[1,c(2,1,6)])
    }
    if (dim(Ki)[1]!=0) {
      print("KiNOTNULL")
      Select_OneMol_Switch=rbind(Select_OneMol_Switch, Ki[1,c(2,1,6)])
    }
  }
}

Select_OneMol_Switch=unique(Select_OneMol_Switch)
Select_OneMol=unique(Select_OneMol)

UI_ToTarget=rbind(Select_OneMol,Select_OneMol_Switch)
unique(UI_ToTarget[,1])
#Linking  mol and scaffolds from  the pareto and the new 100 one

AllInfo5k_Publi=semi_join(Scaff5k_Redone, ParetoSub43,"Molecule")
AllInfo5k_Publi=unique(AllInfo5k_Publi[,c(2,5)])

colnames(AllInfo5k_Publi)[1]="Target"
UI_ToTarget=UI_ToTarget[,-3]
All_Dataset=rbind(AllInfo5k_Publi,UI_ToTarget)
#write.csv(All_Dataset, "AllInfo_Data.csv", row.names = F)

# LOAD CSV WITH HEADERS FROM 'file:///AllInfo_Data.csv' as row FIELDTERMINATOR ','
# MATCH (:ProteinClass)<-[v:memberOf]-(u:UniprotInter)<-[i:memberOf]-(t:Target)<-[a:actifConf9]-(m:Molecule)-[:structureOf]->(c:CompoundName)
# WHERE u.Name=row.Target  AND m.InchiKey = row.Molecule AND c.DBFrom = "CHEMBL"
# WITH u.Name as UniprotInter, t.Name as Target, m.InchiKey as InchiKey, c.MolID as chEMBLID, m.Smiles as SMILES
# RETURN DISTINCT chEMBLID,InchiKey, SMILES, Target
#save as All_Data.csv

#We checked the UI_ToTarget for molecules that would hit the same target for a Human and Mouse/Rat to only keep the Human version.
#We removed them in the All_Data.csv before opening it
Data_Final5k100=read.csv("All_Data.csv", header = T)
Data_Final5k100 <- aggregate(Target~., Data_Final5k100, paste, collapse=";")
#write.csv(Data_Final5k100, "FinalData_5k100.csv", row.names=F)

