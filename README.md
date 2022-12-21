### msPHATE-Seurat

#### What is this?
I integrated multiscale PHATE with Seurat, using scRNA-seq data on the developing adrenal medulla from [Jansky et al.](https://www.nature.com/articles/s41588-021-00806-1), and analysed the results as best I could

#### Orientation
This was my first serious programming project, done as part of a [Vacation Research Scholarship](https://www.unisa.edu.au/research/degrees/scholarships/vacation-research-scholarships/#clinical-health) with UniSA in Adelaide

Please read the final report for context and presentation of the results I decided to keep (including a TLDR/abstract!) [TODO]

#### Scripts
msphate.py is the python script to run multiscale PHATE
msphate.Rmd contains the scripts to export data to multiscale PHATE, import its outputs, and visualise them with Seurat

#### Misc
reports contains the draft analyses I ran in the project (note, it's the woods out there, things may not always work and it's only provided for reference)
It also contains some of the plots generated and necessary to understand what the code is doing (I didn't knit them at the time of making) [TODO]
Chronologically, the "reports" are msphate_report, report_2kgenes, kat_report, final_report, final_report_plots, presentation_plots, combined_report, combined_2
If you are so inclined, here is the Drive folder I was using for notes in the first two weeks of the project

#### Data access
The Seurat objects containing results of all analyses can be found on the University of Adelaide box [TODO]