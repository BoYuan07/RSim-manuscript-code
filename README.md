# RSim-manuscript-code
This is the repository archiving code and data for RSim manuscript--RSim: A Reference-Based Normalization Method via Rank Similarity.

# File Introduction
The "data" folder contains three data sets we use in the manuscript.

The "code" folder contains the RScript and Rmarkdown files for the algorithm and results in the manuscript. Algorithm.R includes the code for RSim and other normalization methods. Files for simulations based on Yan's data can be found in folders named by corresponding downstream analysis. RMarkdown file for real data analysis based on Vangay's data set and Caporaso's data set can be found in folder 'RealData'. All the simulation results can be found in folder 'table'. Figures can be found in folder 'figure'.

# Results Replication
First run code in Algorithm.R, then results can be replicated by corresponding files.

# Session Info
R version 4.2.2 (2022-10-31)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.5

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] compiler_4.2.2  fastmap_1.1.0   cli_3.6.0       htmltools_0.5.4
 [5] tools_4.2.2     rstudioapi_0.14 yaml_2.3.6      rmarkdown_2.19 
 [9] knitr_1.41      xfun_0.36       digest_0.6.31   rlang_1.0.6    
[13] evaluate_0.19  
