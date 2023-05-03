### msPHATE-Seurat

#### What is this?
I integrated **multiscale PHATE** with **Seurat**, using scRNA-seq data on the developing adrenal medulla from [Jansky et al.](https://www.nature.com/articles/s41588-021-00806-1), and analysed the results as best I could. Please see the final report [here](https://htmlpreview.github.io/?https://github.com/artemiyhussnain/msPHATE-Seurat/blob/main/FinalReport_final.html) (this is an html preview). 

This was my first serious programming project, done as part of a [Vacation Research Scholarship](https://www.unisa.edu.au/research/degrees/scholarships/vacation-research-scholarships/#clinical-health) with UniSA in Adelaide

Please read the **final report** for context and presentation of the results I decided to keep

**Scripts** below explains the most important outcome of this project - the three scripts integrating msPHATE into Seurat


#### Scripts

- msphate_pre.R reads in a Seurat object and exports data for use by multiscale PHATE
- msphate.py is the python script to run multiscale PHATE
- msphate_post.R adds the outputs of multiscale PHATE to the Seurat object and visualises the results
- msphate.Rmd contains this basic workflow and scripts for exporting data for all genes, exporting unscaled data, and zooming into clusters


#### Misc

- **reports** contains the draft analyses I ran in the project (note, it's the woods out there, things may not always work and it's only provided for reference)
- The plots generated when making these reports and necessary to understand what the code is doing are all here: https://drive.google.com/drive/folders/1iYqhl6Bo1aZzl_qwWsFtCxqgRU7_FE2I?usp=sharing (I didn't knit them at the time of making...)
- Chronologically, the "reports" are msphate_report, report_2kgenes, kat_report, final_report, final_report_plots, presentation_plots, combined_report, combined_2
- If you are so inclined, here is the Drive folder I was using for notes in the first two weeks of the project: https://drive.google.com/drive/folders/1CYxXftstTIr3SndVfi6KhPeL67zt3oDU?usp=sharing


#### Data access
The Seurat objects containing results of all analyses can be found in https://universityofadelaide.box.com/s/p8aptm3kxksoq4ok9alon4pmm6mj50rl 
