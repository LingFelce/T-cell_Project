##### CCB CLUSTER #####
# for single cell analysis

> sessionInfo()
R version 4.0.1 (2020-06-06)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS:   /Filers/package/R-base/4.0.1/lib64/R/lib/libRblas.so
LAPACK: /Filers/package/R-base/4.0.1/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
 [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] harmony_1.0              Rcpp_1.0.5               plyr_1.8.6               forcats_0.5.0           
 [5] stringr_1.4.0            purrr_0.3.4              readr_1.3.1              tidyr_1.1.2             
 [9] tibble_3.0.3             tidyverse_1.3.0          data.table_1.13.0        ggplot2_3.3.2           
[13] cowplot_1.1.0            patchwork_1.0.0          dplyr_1.0.2              pbmcsca.SeuratData_3.0.0
[17] pbmc3k.SeuratData_3.1.4  panc8.SeuratData_3.0.2   SeuratData_0.2.1         SeuratDisk_0.0.0.9013   
[21] sctransform_0.3.1.9003   Seurat_3.9.9.9010       

loaded via a namespace (and not attached):
  [1] Rtsne_0.15            colorspace_1.4-1      deldir_0.1-25         ellipsis_0.3.1       
  [5] ggridges_0.5.2        fs_1.5.0              rstudioapi_0.11       spatstat.data_1.4-3  
  [9] leiden_0.3.3          listenv_0.8.0         ggrepel_0.8.2         bit64_0.9-7          
 [13] lubridate_1.7.9       fansi_0.4.1           xml2_1.3.2            codetools_0.2-16     
 [17] splines_4.0.1         knitr_1.29            polyclip_1.10-0       jsonlite_1.7.0       
 [21] packrat_0.5.0         broom_0.7.0           ica_1.0-2             dbplyr_1.4.4         
 [25] cluster_2.1.0         png_0.1-7             uwot_0.1.9.9000       shiny_1.5.0          
 [29] compiler_4.0.1        httr_1.4.2            backports_1.1.9       assertthat_0.2.1     
 [33] Matrix_1.2-18         fastmap_1.0.1         lazyeval_0.2.2        cli_2.0.2            
 [37] later_1.1.0.1         htmltools_0.5.0       tools_4.0.1           rsvd_1.0.3           
 [41] igraph_1.2.5          gtable_0.3.0          glue_1.4.2            RANN_2.6.1           
 [45] reshape2_1.4.4        rappdirs_0.3.1        spatstat_1.64-1       cellranger_1.1.0     
 [49] vctrs_0.3.4           nlme_3.1-148          lmtest_0.9-37         xfun_0.16            
 [53] globals_0.13.1        rvest_0.3.6           mime_0.9              miniUI_0.1.1.1       
 [57] lifecycle_0.2.0       irlba_2.3.3           goftest_1.2-2         future_1.20.1        
 [61] MASS_7.3-51.6         zoo_1.8-8             scales_1.1.1          hms_0.5.3            
 [65] promises_1.1.1        spatstat.utils_1.17-0 parallel_4.0.1        RColorBrewer_1.1-2   
 [69] reticulate_1.16       pbapply_1.4-2         gridExtra_2.3         rpart_4.1-15         
 [73] stringi_1.4.6         rlang_0.4.7           pkgconfig_2.0.3       matrixStats_0.56.0   
 [77] lattice_0.20-41       ROCR_1.0-11           tensor_1.5            htmlwidgets_1.5.1    
 [81] bit_1.1-15.2          tidyselect_1.1.0      parallelly_1.21.0     RcppAnnoy_0.0.16     
 [85] magrittr_1.5          R6_2.4.1              generics_0.0.2        DBI_1.1.0            
 [89] haven_2.3.1           pillar_1.4.6          withr_2.2.0           mgcv_1.8-31          
 [93] fitdistrplus_1.1-1    survival_3.1-12       abind_1.4-5           future.apply_1.5.0   
 [97] modelr_0.1.8          crayon_1.3.4          hdf5r_1.3.3           KernSmooth_2.23-17   
[101] plotly_4.9.2.1        readxl_1.3.1          grid_4.0.1            blob_1.2.1           
[105] reprex_0.3.0          digest_0.6.25         xtable_1.8-4          httpuv_1.5.4         
[109] munsell_0.5.0         viridisLite_0.3.0  


# for pseudo/minibulk analysis

> sessionInfo()
R version 4.0.1 (2020-06-06)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS:   /Filers/package/R-base/4.0.1/lib64/R/lib/libRblas.so
LAPACK: /Filers/package/R-base/4.0.1/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
 [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] clusterProfiler_3.16.0      ggrepel_0.8.2               RColorBrewer_1.1-2         
 [4] DT_0.13                     forcats_0.5.0               stringr_1.4.0              
 [7] purrr_0.3.4                 readr_1.3.1                 tidyr_1.1.2                
[10] tibble_3.0.3                tidyverse_1.3.0             data.table_1.13.0          
[13] genefilter_1.70.0           gplots_3.0.3                BiocParallel_1.22.0        
[16] statmod_1.4.34              edgeR_3.30.3                limma_3.44.3               
[19] stringi_1.4.6               ggplot2_3.3.2               dplyr_1.0.2                
[22] pheatmap_1.0.12             DESeq2_1.28.1               SummarizedExperiment_1.18.2
[25] DelayedArray_0.14.1         matrixStats_0.56.0          GenomicRanges_1.40.0       
[28] GenomeInfoDb_1.24.2         org.Hs.eg.db_3.11.4         AnnotationDbi_1.50.1       
[31] IRanges_2.22.2              S4Vectors_0.26.1            Biobase_2.48.0             
[34] BiocGenerics_0.34.0        

loaded via a namespace (and not attached):
  [1] reticulate_1.16        tidyselect_1.1.0       RSQLite_2.2.0          htmlwidgets_1.5.1     
  [5] grid_4.0.1             Rtsne_0.15             scatterpie_0.1.4       munsell_0.5.0         
  [9] codetools_0.2-16       ica_1.0-2              future_1.20.1          miniUI_0.1.1.1        
 [13] withr_2.2.0            GOSemSim_2.14.0        colorspace_1.4-1       knitr_1.29            
 [17] rstudioapi_0.11        Seurat_3.9.9.9010      ROCR_1.0-11            tensor_1.5            
 [21] DOSE_3.14.0            listenv_0.8.0          labeling_0.3           rstan_2.21.2          
 [25] urltools_1.7.3         GenomeInfoDbData_1.2.3 polyclip_1.10-0        farver_2.0.3          
 [29] bit64_0.9-7            downloader_0.4         parallelly_1.21.0      vctrs_0.3.4           
 [33] generics_0.0.2         xfun_0.16              R6_2.4.1               graphlayouts_0.7.0    
 [37] rsvd_1.0.3             locfit_1.5-9.4         gridGraphics_0.5-0     fgsea_1.14.0          
 [41] bitops_1.0-6           spatstat.utils_1.17-0  assertthat_0.2.1       promises_1.1.1        
 [45] scales_1.1.1           ggraph_2.0.3           enrichplot_1.8.1       gtable_0.3.0          
 [49] globals_0.13.1         processx_3.4.3         goftest_1.2-2          tidygraph_1.2.0       
 [53] rlang_0.4.7            splines_4.0.1          lazyeval_0.2.2         europepmc_0.4         
 [57] broom_0.7.0            inline_0.3.15          BiocManager_1.30.10    yaml_2.2.1            
 [61] reshape2_1.4.4         abind_1.4-5            modelr_0.1.8           crosstalk_1.1.0.1     
 [65] backports_1.1.9        httpuv_1.5.4           qvalue_2.20.0          tools_4.0.1           
 [69] ggplotify_0.0.5        ellipsis_0.3.1         ggridges_0.5.2         Rcpp_1.0.5            
 [73] plyr_1.8.6             progress_1.2.2         zlibbioc_1.34.0        RCurl_1.98-1.2        
 [77] ps_1.3.4               prettyunits_1.1.1      rpart_4.1-15           deldir_0.1-25         
 [81] viridis_0.5.1          pbapply_1.4-2          cowplot_1.1.0          zoo_1.8-8             
 [85] haven_2.3.1            cluster_2.1.0          fs_1.5.0               magrittr_1.5          
 [89] DO.db_2.9              triebeard_0.3.0        lmtest_0.9-37          reprex_0.3.0          
 [93] RANN_2.6.1             packrat_0.5.0          fitdistrplus_1.1-1     evaluate_0.14         
 [97] hms_0.5.3              patchwork_1.0.0        mime_0.9               xtable_1.8-4          
[101] XML_3.99-0.3           readxl_1.3.1           gridExtra_2.3          compiler_4.0.1        
[105] KernSmooth_2.23-17     V8_3.2.0               crayon_1.3.4           StanHeaders_2.21.0-5  
[109] htmltools_0.5.0        mgcv_1.8-31            later_1.1.0.1          geneplotter_1.66.0    
[113] RcppParallel_5.0.1     lubridate_1.7.9        DBI_1.1.0              tweenr_1.0.1          
[117] dbplyr_1.4.4           MASS_7.3-51.6          rappdirs_0.3.1         Matrix_1.2-18         
[121] cli_2.0.2              gdata_2.18.0           igraph_1.2.5           pkgconfig_2.0.3       
[125] rvcheck_0.1.8          plotly_4.9.2.1         xml2_1.3.2             annotate_1.66.0       
[129] XVector_0.28.0         rvest_0.3.6            callr_3.4.3            digest_0.6.25         
[133] sctransform_0.3.1.9003 RcppAnnoy_0.0.16       spatstat.data_1.4-3    fastmatch_1.1-0       
[137] rmarkdown_2.3          cellranger_1.1.0       leiden_0.3.3           uwot_0.1.9.9000       
[141] curl_4.3               shiny_1.5.0            gtools_3.8.2           lifecycle_0.2.0       
[145] nlme_3.1-148           jsonlite_1.7.0         viridisLite_0.3.0      fansi_0.4.1           
[149] pillar_1.4.6           lattice_0.20-41        loo_2.3.1              GO.db_3.11.4          
[153] fastmap_1.0.1          httr_1.4.2             pkgbuild_1.1.0         survival_3.1-12       
[157] glue_1.4.2             spatstat_1.64-1        png_0.1-7              bit_1.1-15.2          
[161] ggforce_0.3.1          blob_1.2.1             caTools_1.18.0         memoise_1.1.0         
[165] irlba_2.3.3            future.apply_1.5.0  



##### BMRC CLUSTER #####

> sessionInfo()
R version 3.6.2 (2019-12-12)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /gpfs3/apps/eb/skylake/software/OpenBLAS/0.3.7-GCC-8.3.0/lib/libopenblas_skylakexp-r0.3.7.so

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] harmony_1.0          Rcpp_1.0.5           forcats_0.4.0        stringr_1.4.0        purrr_0.3.4         
 [6] readr_1.3.1          tidyr_1.1.2          tibble_3.0.4         tidyverse_1.3.0      data.table_1.12.8   
[11] ggplot2_3.3.3        cowplot_1.1.1        patchwork_1.1.0.9000 dplyr_1.0.2          plyr_1.8.5          
[16] Seurat_3.2.3        

loaded via a namespace (and not attached):
  [1] Rtsne_0.15            colorspace_2.0-0      deldir_0.1-23         ellipsis_0.3.1        ggridges_0.5.1       
  [6] fs_1.3.1              rstudioapi_0.13       spatstat.data_1.4-0   leiden_0.3.6          listenv_0.8.0        
 [11] npsurv_0.4-0          ggrepel_0.8.1         fansi_0.4.1           lubridate_1.7.4       xml2_1.2.2           
 [16] codetools_0.2-16      splines_3.6.2         lsei_1.2-0            knitr_1.26            polyclip_1.10-0      
 [21] jsonlite_1.7.2        broom_0.5.3           ica_1.0-2             cluster_2.1.0         dbplyr_1.4.2         
 [26] png_0.1-7             uwot_0.1.10           shiny_1.4.0           sctransform_0.3.2     compiler_3.6.2       
 [31] httr_1.4.1            backports_1.1.5       assertthat_0.2.1      Matrix_1.2-18         fastmap_1.0.1        
 [36] lazyeval_0.2.2        cli_2.2.0             later_1.0.0           htmltools_0.4.0       tools_3.6.2          
 [41] rsvd_1.0.3            igraph_1.2.4.2        gtable_0.3.0          glue_1.4.2            RANN_2.6.1           
 [46] reshape2_1.4.3        spatstat_1.62-2       scattermore_0.7       cellranger_1.1.0      vctrs_0.3.6          
 [51] gdata_2.18.0          nlme_3.1-143          lmtest_0.9-37         xfun_0.11             globals_0.14.0       
 [56] rvest_0.3.5           mime_0.7              miniUI_0.1.1.1        lifecycle_0.2.0       irlba_2.3.3          
 [61] gtools_3.8.1          goftest_1.2-2         future_1.21.0         MASS_7.3-51.4         zoo_1.8-6            
 [66] scales_1.1.1          hms_0.5.2             promises_1.1.0        spatstat.utils_1.15-0 parallel_3.6.2       
 [71] RColorBrewer_1.1-2    reticulate_1.13       pbapply_1.4-2         gridExtra_2.3         rpart_4.1-15         
 [76] stringi_1.4.3         caTools_1.17.1.3      rlang_0.4.10          pkgconfig_2.0.3       matrixStats_0.57.0   
 [81] bitops_1.0-6          lattice_0.20-38       ROCR_1.0-7            tensor_1.5            htmlwidgets_1.5.1    
 [86] tidyselect_1.1.0      parallelly_1.23.0     RcppAnnoy_0.0.18      magrittr_2.0.1        R6_2.5.0             
 [91] gplots_3.0.1.1        generics_0.1.0        DBI_1.1.0             pillar_1.4.7          haven_2.2.0          
 [96] withr_2.3.0           mgcv_1.8-31           fitdistrplus_1.0-14   survival_3.1-8        abind_1.4-5          
[101] future.apply_1.7.0    modelr_0.1.5          crayon_1.3.4          KernSmooth_2.23-16    plotly_4.9.1         
[106] readxl_1.3.1          grid_3.6.2            reprex_0.3.0          digest_0.6.27         xtable_1.8-4         
[111] httpuv_1.5.2          munsell_0.5.0         viridisLite_0.3.0    
