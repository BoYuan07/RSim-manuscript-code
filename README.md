# RSim-manuscript-code
This is the repository archiving code and data for RSim manuscript--RSim: A Reference-Based Normalization Method via Rank Similarity.

# File Introduction
The "data" folder contains three data sets we use in the manuscript.

The "code" folder contains the RScript and Rmarkdown files for the algorithm and results in the manuscript. Algorithm.R includes the code for RSim and other normalization methods. Files for simulations based on Yan's data can be found in folders named by corresponding downstream analysis. RMarkdown file for real data analysis based on Vangay's data set and Caporaso's data set can be found in folder 'RealData'. All the simulation results can be found in folder 'table'. Figures can be found in folder 'figure'.

# Results Replication
First run code in Algorithm.R, then results can be replicated by corresponding files.

# Session Info
R version 4.1.1 (2021-08-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux Server 7.9 (Maipo)

Matrix products: default
BLAS/LAPACK: /usr/local/src/openblas/0.3.12/gcc/Sandy.Bridge/lib/libopenblas_sandybridgep-r0.3.12.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.1.1
