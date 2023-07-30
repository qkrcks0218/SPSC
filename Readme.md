# Replication Files for Park and Tchetgen Tchetgen (2023)

This Github repository contains replication files for ``Single Proxy Synthetic Control'' [Park and Tchetgen Tchetgen (2023)](https:// "SPSC") (SPSC).


## Folders

### 1.Simulation

1.Simulation folder contains 0.Function_SSC_Simulation.R, 1.Simulation.R, and 2.PaperReport.R files. These files replicate the simulation study in Section 4 of the paper. 

* 0.Function_SSC_Simulation.R contains functions used for the simulation study.

* 1.Simulation.R is used to replicate the simulation study with a parameter BATCH$`\in \{1,\ldots,6000\}`$. For each BATCH parameter, $`d \in \{0,1\}`$ (0 = without covariate, 1 = with covariate), $`b \in \{1,2\}`$ (1 = constant average treatment effect on the treated (ATT), 2 = linear ATT), $`n \in \{02,05,09\}`$ (number of donors), $`i \in \{001,\ldots,500\}`$ (simulation repetition index) is assigned, and 4 csv files with a format SSC\_T$`[t]`$\_D$`[d]`$\_B$`[b]`$\_N $`[n]`$\_Iter$`[i]`$.csv ($`t \in \{00050,00100,00250,01000\}`$) are generated in Result\_D$`[d]`$_B$`[b]`$\_N$`[n]`$ folder. It is recommended to use a parallel computing system by parallelizing the BATCH parameter. 

* 2.PaperReport.R is used to summarize the results into tables and figures of the paper and supplementary material.

### 1.Simulation_Supp

1.Simulation_Supp folder contains 0.Function_SSC_Simulation.R, 1.Simulation.R, and 2.PaperReport.R files. These files replicate the simulation study in Section 4 of the paper. 

* 0.Function_SSC_Simulation2.R contains functions used for the simulation study.

* 1.NoAdjustment.R and 2.NoAdjustment.R are used to replicate the simulation study in Section S1.9 of the supplementary material based on the simulation scenario given in [Cattaneo et al. (2021)](https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1979561, "SCPI"). The two code files have similar structures except that 1.NoAdjustment.R uses time-invariant $`g`$ and 2.NoAdjustment.R uses time-varying $`g_t`$. These R codes involve with a parameter $`i \in \{1,\ldots,500\}`$, and for each BATCH, a csv file with a format Result_BATCH$`[i]`$.csv is created in Result_Raw or Result_Raw_Time folder. 

* 3.Summary.R is used to summarize the results into tables of the supplementary material.


### 2.Data

The dataset used in the paper is first analyzed in [Fohlin and Lu (2021)](https://www.aeaweb.org/articles?id=10.1257/pandp.20211097, "data"), and is publicly available at the following [[link](https://www.aeaweb.org/articles?id=10.1257/pandp.20211097)].

Data folder contains the following files, which replicate the data analysis results in Section 5 of the paper.

* 0.Function_SPSC_Data.R contains functions used for the data analysis.

* 0.DataCleaning.R and 0.DataCleaning_Placebo.R are used to clean the dataset based on the actual treatment time and the placebo treatment time, respectively. 

* 1-1.ATT_SPSC.R and 1-2.ATT_Time_SPSC.R are used to estimate the ATT based on the SPSC approach with time-invariant $`g`$ and time-varying $`g_t`$, respectively.

* 1-3.ATT_Abadie.R is used to estimate the ATT based on the approach of [Abadie et al. (2010)](https://www.tandfonline.com/doi/abs/10.1198/jasa.2009.ap08746 "Abadie").

* 1-4.ATT_Placebo_SPSC.R is used to estimate the placebo ATT based on the SPSC approach.

* 2-1.Conformal_SPSC.R and 2-2.Conformal_Time_SPSC.R are used to estimate prediction intervals of the ATT based on the conformal inference approach of [Chernozhukov et al. (2021)](https://www.tandfonline.com/doi/abs/10.1080/01621459.2021.1920957, "Conformal") adapted to the SPSC setting.

* 2-3.Conformal_SCPI.R is used to estimate prediction intervals of the ATT based on the approach of [Cattaneo et al. (2021)](https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1979561, "SCPI")

* 2-4.Conformal_Placebo_SPSC.R is used to estimate prediction intervals of the placebo ATT based on the conformal inference approach of [Chernozhukov et al. (2021)](https://www.tandfonline.com/doi/abs/10.1080/01621459.2021.1920957, "Conformal") adapted to the SPSC setting.

* 3-1.Table_ATT.R, 3-2.Plot_Conformal.R, and 3-3.Table_Plot_Placebo.R are used to summarize the results into tables and figures of the paper and supplementary material.

* 10-1.Analysis_OverlapSelection.R, 10-2.Analysis_LassoSelection.R, and 10-3.Analysis_RegressionSelection.R provide details on how donors are selected. 



## References

Alberto Abadie, Alexis Diamond & Jens Hainmueller (2010) **Synthetic Control Methods for Comparative Case Studies: Estimating the Effect of California’s Tobacco Control Program**, _Journal of the American Statistical Association_, 105:490, 493-505 [[link](https://www.tandfonline.com/doi/abs/10.1198/jasa.2009.ap08746)]

Victor Chernozhukov, Kaspar Wüthrich & Yinchu Zhu (2021) **An Exact and Robust Conformal Inference Method for Counterfactual and Synthetic Controls**, _Journal of the American Statistical Association_, 116:536, 1849-1864 [[link](https://www.tandfonline.com/doi/abs/10.1080/01621459.2021.1920957)]

Matias D. Cattaneo, Yingjie Feng & Rocio Titiunik (2021) **Prediction Intervals for Synthetic Control Methods**, _Journal of the American Statistical Association_, 116:536, 1865-1880 [[link](https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1979561)]

Caroline Fohlin & Zhikun Lu (2021) **How Contagious Was the Panic of 1907? New Evidence from Trust Company Stocks** _AEA Papers and Proceedings_, 111: 514-19 [[link](https://www.aeaweb.org/articles?id=10.1257/pandp.20211097)]

Chan Park & Eric Tchetgen Tchetgen (2023) **Single Proxy Synthetic Control**, _arXiv:_ [[link](https:// "SPSC")]