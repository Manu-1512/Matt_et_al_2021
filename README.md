# Matt_et_al_2021
This is a simultaneous analysis of genes and Transposable elements (TEs) from 10X data, using the lovain clustering. We performed the analysis using “Alevin” as a tool integrated with the salmon software that uses EM algorithms for the analysis of 3’ tagged-end single-cell sequencing data.

_**1.**_ The code named ```Preparing_files_for_Alevin_bash.sh``` is step-wise codes where we indexed the concatenated genes and TEs transcriptome and genome reference file using salmon.

_**2.**_ we run Alevin on all the 10X generated samples at once using ```Running_Alevin.sh```

_**3.**_ Next file ```Imorting_Merging_files_fromAlevin_output.r``` is where we import and merge all datasets in one seurat object. We also calculate the percentage of mitochondrial or ribosomal genes, and the statistics of numbers of genes, UMIs in single cells

_**4.**_ Now we remove the doublets, and annotate the cells from merged but not integrated objects using the codes in                           ```Doublet_removed_Annotating_cells_from_merged_object.r```

**DATA ANALYSIS** - 
                Once the datasets were merged or even as individual files, we performed the integration, normalization, dimension reduction, visualization, and                     trajectory analysis of the single cellular clusters. Also, TE loci graphs were generated using the normalized data. All of the relevant codes are                   present in the file    ```Data_integration_Dimension_reduction_visualization.r```

  _**TE family level analysis**_ 
               This is followup of previous code where we analyze the data for TEs at family level. Also by removing all the TEs and comparing with genes only                      analysis. Code is given in ```TE_analysis_at_family_level_and_only_genes_analysis.r```
               
              
