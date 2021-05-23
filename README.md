# Matt_et_al_2021
This is a simultaneous analysis of genes and Transposable elements (TEs) from 10X data, using the lovain clustering. We performed the analysis using “Alevin” as a tool integrated with the salmon software that uses EM algorithms for the analysis of 3’ tagged-end single-cell sequencing data.

First code named ```Preparing_files_for_Alevin_bash.sh``` is step-wise codes where we indexed the concatenated genes and TEs transcriptome and genome reference file using salmon.

Then we run Alevin on all the 10X generated samples at once using ```Running_Alevin.sh```

Next file ```Imorting_Merging_files_fromAlevin_output.r``` is where we import and merge all datasets in one seurat object. We also calculate the percentage of mitochondrial or ribosomal genes, and the statistics of numbers of genes, UMIs in single cells

