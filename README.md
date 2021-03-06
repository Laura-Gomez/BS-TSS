# BS-TSS

This program finds BS-TSS interactions that are in a regulatory neighborhood, and that could have a biological function associated.

## DESCRIPTION

### INPUT

For each input file, a list of columns is provided. All columns must exist. When the information is no available NA is allowed. The mandatory fields are shown in bold letters

**Binding sites from RegulonDB** [regulon.bs]
 - TF ID
 - **TF name**
 - TFBS ID
 - TFBS left position
 - TFBS rigth position 
 - TFBS strand
 - TFBS associated gene ID
 - TFBS associated TU 
 - **TFBS effect**
 - TFBS associated promoter
 - **TFBS distance to TSS**
 - TFBS sequence
 - TFBS evidence
 - **TFBS evidence level**

**Expression data** [expression.data]
 - Run
 - Sample
 - **TF**
 - **Target**
 - **LogFoldFPKM**
 - LogFoldTPM
 - FPKM
 - Counts
 - WildTypeFPKM
 - WildTypeTPM
 - WTCounts

**Peaks data** [peaks.data]
 - Exp
 - sample
 - TF
 - type
 - ID1
 - Start
 - Stop
 - **PeakPos**
 - Height
 - No1
 - No2
 - Shift
 - Bnumber
 - **Gene**
 - Dist

**TSS information from RegulonDB (known TSSs - promoters download file)** [regulon.tss]
 - ID
 - Promoter.Name
 - **TSS**
 - Sigma
 - Strand
 - GI
 - **Gene**
 - PosLeft
 - PosRigth
 - Evidence
 - Bnumber
 - PromSeq
 - DistPromGene

**HT TSSs** [ht.tss]
This code was written using Gisella Storz data
 - **TSSPosition**
 - RPKM
 - Promoter
 - Strand
 - RelPos
 - **Gene**
 - Bnumber
 - LeftGene
 - RigthGene
 - **Orientation**
 - TSSClass
 - Enrichment
 - evidence


## OUTPUT

Each BS is associated with one gene, and each gene is assigned to the category of activated or repressed

### OUTPUT FILES 

 - Histograms of the distance between the known BSs for a specific TF and the known TSS´s associated to those site
 - *TSS.activated.all*: Number of possible TSS-BS interactions for all BSs associated to activated genes. For regulon TSSs
 - *TSS.activated.pass*: Number of TSS-BS interactions that pass the distance filter for all BSs associated to activated genes. For Regulon TSSs
 - *TSS.repressed.all*: Same as TSS.activated.all. BSs associated to repressed genes
 - *TSS.repressed.pass*: Same as TSS.activated.pass. BSs associated to repressed genes
 - *HT-TSS.activated.all*: Same as TSS.activated.all. TSSs taken from HT experiment.
 - *HT-TSS.activated.pass*: Same as TSS.activated.pass. TSSs taken from HT experiment.
 - *HT-TSS.repressed.all*: Same as TSS.repressed.all. TSSs taken from HT experiment.
 - *HT-TSS.repressed.pass*: Same as TSS.repressed.pass. TSSs taken from HT experiment.
	
### OUTPUT DIRECTORIES
Each output directory containd two directories: Activated and Repressed:

 - *Activated*: Contains all TSS-BS distances for all genes classifies as Activated in the expression data
 - *Repressed*:  Contains all TSS-BS distances for all genes classifies as Repressed in the expression data
 
Each subdirectory contains one file per gene. In these files, each columns represents a BS and each row represents a TSS.

The output directories are:
 - *TSS*: All distances from BSs to known TSSs (regulon TSSs)
 - *TSS-PASS*: Distances that passed the distance filter (regulon TSSs)
 - *HT-TSS*: All distances from BSs to HT TSSs (eg, Storsz TSSs)
 - *HT-TSS-PASS*: Distances that passed the distance filter (eg, Storsz TSSs)
 
 
### RUN

 Rscript --vanilla Distance_TSS_BS.R regulon.bs TF expression.data peaks.data regulon.tss ht.tss Example-Data/INPUT
 
