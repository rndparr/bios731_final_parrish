# bios731_final_parrish


---

- [Directories \& Files](#directories-and-files)
- [Reproducibility](#reproducibility)
- [R Session Info](#R-session-info)
- [LaTeX Info](#latex-info)


---


## Directories and Files


```
.
├── README.md
├── data
│   ├── ProteinAbundance_ROSMAP_annot.txt
│   ├── TranscriptAbundance_ROSMAP_annot.txt
│   ├── cis_PWAS_results.txt
│   ├── cis_TWAS_results.txt
│   ├── cistrans_PWAS_results.txt
│   ├── cistrans_TWAS_results.txt
│   ├── cisPWAS_pathDIP_Enrichment_analysis_for_Panther_Pathway.txt
│   ├── cisTWAS_pathDIP_Enrichment_analysis_for_Panther_Pathway.txt
│   ├── cistransPWAS_pathDIP_Enrichment_analysis_for_Panther_Pathway.txt
│   └── cistransTWAS_pathDIP_Enrichment_analysis_for_Panther_Pathway.txt
├── drawio_figs
│   ├── bio_bg_2-1.drawio.pdf
│   ├── bio_bg_2-GWAS_2bg.pdf
│   └── bio_bg_2-Two_Step_xWAS_1.drawio.pdf
├── figs
│   ├── manplots_report.png
│   ├── manplots_slides.png
│   ├── paths_plot.pdf
│   ├── paths_venn_report.pdf
│   ├── venn_report.pdf
│   ├── venn_slides.pdf
│   └── xWAS.jpg
├── source
│   ├── asso_study_results.R
│   ├── manhattan_plot.R
│   └── pathDIP_results.R
├── tables
│   ├── cis-PWAS_sig_genes.tex
│   ├── cis-TWAS_sig_genes.tex
│   ├── cistrans-PWAS_sig_genes.tex
│   ├── cistrans-TWAS_sig_genes.tex
│   ├── overlap_paths.tex
│   ├── summary_report.tex
│   └── summary_slides.tex
└── tex
    ├── abbrvnat-maxbibnames4_sortnumeric.bst
    ├── biblio.bib
    ├── report.pdf
    ├── report.tex
    ├── slides_content.tex
    ├── slides_formatting.tex
    ├── slides_handout.pdf
    ├── slides_handout.tex
    ├── slides_present.pdf
    └── slides_present.tex
```

- `tex/` contains the report and slide tex files `report.tex`, `slides_present.tex`, and `slides_handout.tex`, as well as the .pdf output of compiling these files, the bibliography information file `biblio.bib`, the bibilography format specification file `abbrvnat-maxbibnames4_sortnumeric.bst`, slide formatting file `slides_formatting.tex`, and slides content file `slides_content.tex`
- `data/` contains the combined output of all genes for each of the four xWASs in the `*_results.txt` files, the output of the pathDIP enrichment analysis files `*_pathDIP_Enrichment_analysis_for_Panther_Pathway.txt`, and two annotation files `ProteinAbundance_ROSMAP_annot.txt` and `TranscriptAbundance_ROSMAP_annot.txt`.
- `source/` contains R-scripts used for this report: `asso_study_results.R`, `pathDIP_results.R`, and `manhattan_plot.R`.
- `tables/` contains tex files of tables output by the `asso_study_results.R` and `pathDIP_results.R` scripts
- `figs/` contains figure files output by the `asso_study_results.R` and `pathDIP_results.R` scripts as well as some additional figures
- `drawio_figs/` contains additional chart figures


---


## Reproducibility


### Required R packages (local)

- cowplot
- ggnewscale
- ggplot2
- ggrepel
- ggVennDiagram
- here
- reshape2
- reshape2
- stringr
- xtable

```R
install.packages(c("cowplot", "ggnewscale", "ggplot2", "ggrepel", "ggVennDiagram", "here", "reshape2", "reshape2", "stringr", "xtable"))
```

### Making Tables and Figures

```bash
Rscript ./source/asso_study_results.R
Rscript ./source/pathDIP_results.R
# note that manhattan_plots.R is used by asso_study_results.R but should not be called directly
```

### Compiling the Report and Slides

The report and slides must be compiled with LaTeX, biber run on this output, and then run LaTeX again to compile a PDF with bibliography.

```bash
cd tex

pdflatex -interaction=nonstopmode report.tex
pdflatex -interaction=nonstopmode slides_handout.tex
pdflatex -interaction=nonstopmode slides_present.tex

biber report
biber slides_handout
biber slides_present

pdflatex -interaction=nonstopmode report.tex
pdflatex -interaction=nonstopmode slides_handout.tex
pdflatex -interaction=nonstopmode slides_present.tex
```


---

## R Session Info

```R
R version 4.5.2 (2025-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] xtable_1.8-8          reshape2_1.4.5        here_1.0.2           
[4] ggVennDiagram_1.5.7   ggnewscale_0.5.2.9000 ggrepel_0.9.8        
[7] ggplot2_4.0.2        

loaded via a namespace (and not attached):
 [1] gtable_0.3.6       dplyr_1.2.1        compiler_4.5.2     tidyselect_1.2.1  
 [5] Rcpp_1.1.1         stringr_1.6.0      systemfonts_1.3.2  scales_1.4.0      
 [9] textshaping_1.0.5  R6_2.6.1           plyr_1.8.9         labeling_0.4.3    
[13] generics_0.1.4     tibble_3.3.1       rprojroot_2.1.1    pillar_1.11.1     
[17] RColorBrewer_1.1-3 rlang_1.2.0        stringi_1.8.7      S7_0.2.1          
[21] cli_3.6.6          withr_3.0.2        magrittr_2.0.5     grid_4.5.2        
[25] cowplot_1.2.0      lifecycle_1.0.5    vctrs_0.7.3        glue_1.8.1        
[29] farver_2.1.2       ragg_1.5.2         tools_4.5.2        pkgconfig_2.0.3   
```

---

## LaTeX Info

```
biber version: 2.21

pdfTeX 3.141592653-2.6-1.40.29 (TeX Live 2026)
kpathsea version 6.4.2
Copyright 2026 Han The Thanh (pdfTeX) et al.
There is NO warranty.  Redistribution of this software is
covered by the terms of both the pdfTeX copyright and
the Lesser GNU General Public License.
For more information about these matters, see the file
named COPYING and the pdfTeX source.
Primary author of pdfTeX: Han The Thanh (pdfTeX) et al.
Compiled with libpng 1.6.55; using libpng 1.6.55
Compiled with zlib 1.3.2; using zlib 1.3.2
Compiled with xpdf version 4.06
```
