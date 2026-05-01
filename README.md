# bios731_final_parrish


---

- [Directories \& Files](#directories-and-files)
- [Reproducibility](#reproducibility)
- [LaTeX Info](#latex-info)


---


## Directories and Files


```
.
├── README.md
├── drawio_figs
│    ├── bio_bg_2-1.drawio.pdf
│    ├── bio_bg_2-GWAS_2bg.pdf
│    ├── bio_bg_2-Two_Step_xWAS_1.drawio.pdf
├── figs
│    ├── manplots_report.png
│    ├── manplots_slides.png
│    ├── paths_plot.png
│    ├── paths_venn_report-crop.pdf
│    ├── venn-crop.pdf
│    ├── venn_report-crop.pdf
│    └── xWAS.jpg
├── final.sublime-project
├── final.sublime-workspace
└── tex
    ├── abbrvnat-maxbibnames4_sortnumeric.bst
    ├── biblio.bib
    ├── report.pdf
    ├── report.tex
    ├── report_figures.tex
    ├── report_tables.tex
    ├── slides_content.tex
    ├── slides_formatting.tex
    ├── slides_handout.pdf
    ├── slides_handout.tex
    ├── slides_present.pdf
    └── slides_present.tex
```

---


## Reproducibility

The report and slides must be compiled with LaTeX, biber run on this output, and then run LaTeX again to compile a PDF with bibliography.

```bash
latex report.tex
latex slides_content.tex

biber report
biber slides_content

latex report.tex
latex slides_content.tex
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

