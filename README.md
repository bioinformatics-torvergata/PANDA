# PAN-cancer Data Analysis Web Tool (PANDA)
## Bioinformatics tools for statistical analysis of tumor gene expression

PANDA (PAN-cancer Data Analysis web tool) is a bioinformatics web server that provides tools for analyzing RNAseq data collected by TCGA. The available analyses leverage clinical data to perform differential gene expression and overall survival analyses with different options.


https://panda.bio.uniroma2.it


## Analyses available:

### **Transcriptomic Analyses**
Perform differential expression and correlation analyses by selecting tumor type, clinical features, gene, and miRNA:
- **DESeq2 Analysis**: Differential Expression Analysis for a single tumor.
- **Pan-Cancer Differential Expression Analysis**.
- **Boxplot Protein**: Visualization of expression data.

### **Proteomic Analyses**
Differential expression analysis based on tumor type, clinical features, and protein name:
- **Differential Expression Analysis Single Tumor**.
- **Pan-Cancer Differential Expression Analysis**.

### **Mutational Analyses**
Perform mutational data analysis by selecting tumor type or gene of interest:
- **Tumor Mutation Analysis**.
- **Oncoplot**.
- **Somatic Interaction Analysis**.
- **Gene Mutation Analysis**.
- **Differential Expression for Mutated Status (DESeq2)**.
- **Differentially Mutated Gene by Clinical Feature**.

### **Survival Analysis**
Survival analysis using mutational, transcriptomic, and proteomic data based on tumor and clinical characteristics:
- **Overall Survival**.
- **Overall Survival with Pathway Activity Score**.
- **Overall Survival with Gene Mutation Status**.

### **Cell Types Analyses**
Cell mixture deconvolution estimates the proportions of different cell types in a mixed sample using gene expression data:
- **Cell-mixture Deconvolution**.
- **Correlation between Cell Type and Pathways**.

## How to Cite PANDA
If you use PANDA, even just for preliminary analysis, please cite: 


## Contact
For questions, support, bug reports, proposals, cooperations or any other communication, feel free to contact:

**Gerardo Pepe** (Conceptualized and designed the study, supervised the development of the webserver and implemented the scripts for data analysis tools): gerardo.pepe@uniroma2.it

**Chiara Notturno Granieri** (Implemented both back-end and front-end components and scripts for data analysis tools): chiara.notturno@gmail.com

**Romina Appierdo** (Worked on the implementation of scripts for data analysis tools): romina.appierdo@uniroma2.it

