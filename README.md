# TCR_IMGT

Jupyter+R [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mattobu83/TCR_IMGT/master?filepath=IMGT_excel_processing.ipynb)

To run the R script, run the following in the terminal: 

```bash
Rscript TCR_IMGT_excel.R In.xls Out.txt In2.xls
```

where In.xls is the excel output from VQuest and Out.txt is the preferred output file name. A second input file can be supplied,In2.xls,if there were more than 50 sequences in the run (Vquest can only process 50 sequences at a time). If there is no second file, the third argument can be omitted. 
