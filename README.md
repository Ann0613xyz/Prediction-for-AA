Please uniformly change the path "C:\Users\anxin\Desktop" in all codes to the path where your current “integrated_model” folder is located.

If you want to reproduce the results, it is recommended to run the various code files in the following order:

1.annotation.R in Folder GSE68801, GSE148346,…【This step is for initial data processing. Before running this step, data acquisition is required: create a folder named "GSE148346" in the current directory and open it. GPL570-55999.txt and GSE148346_series_matrix.txt.gz are downloaded from the corresponding data links on GEO. After decompressing GSE148346_series_matrix.txt.gz, the information is manually split into two files: clinical information and GSE148346_series_matrix.txt. GSE68801 and GSE148346 are uploaded as example datasets, and the same applies to the other datasets.】

2.combine.R：Merge datasets and record merging results (generating “all.txt”, “diff.txt”, “normalize.txt”, “diffGeneExp.txt”), PCA plots and box plots (Supplementary Figure 1), DEG analysis and heatmap for train set (Figure 1b) (allcontrol and allAA mentioned in the code are the sample names of the Con group and AA group respectively, obtained from clinical information).

3.test_process.R：Generating “testall.txt”, “testdiff.txt”, “testnormalize.txt”, “testdiffGeneExp.txt”), DEG analysis and heatmap for test set (Figure 1d).

4.volcano_plot.R: Volcano plots in train set and test set respectively (Figure 1a,c).

5.venn of DEG.R: Intersected DEGs between train and test set (Figure 1e,f).

6.GO.R: GO analysis for DEGs in train set (Figure 2c,d).

7.GOtest.R: GO analysis for DEGs in test set (Figure 2e,f).

8.KEGG.R: KEGG analysis for DEGs in train set (Figure 2a,b).

9.GSEA.R: GSEA for all genes in train set (Figure 3a,b).

10.GSEAtest.R: GSEA for all genes in test set (Figure 3c,d)

11.CIBERSORTrun.R: Generating “CIBERSORT-Results.txt” (Result of immune filtrated in train set).

12.CIBERSORTrun-test.R: Generating “CIBERSORT- testResults.txt” (Result of immune filtrated in test set).

13.cibersort_immune_barplot.R: Generating proportion of immune cells for train set, correlation between immune cells in train set, correlation between immune cells in test set

14.cibersort_vioplot.R: proportion of immune cells’ difference between Con and AA group in train set. (Figure 3h)

15,lasso.R: Lasso Regression for feature genes selection and bootstrap resampling validation (Figure 4 a,b,j)

16.randomforest.R: RF for feature genes selection and bootstrap resampling validation (Figure 4 c,d,k)

17.venn and ROC.R: ROC curves for three feature genes (Figure 4 f,g,h)

18.feature_genes_vioplot.R: Vioplots demonstrating down-regulation of feature genes (Figure 4i)

19.cibersort_immuneCor.R: Demonstrating correlation between immune genes and feature genes. (Figure 4o,p,q)

20.cibersort_corimmune_Lollipop.R: Visualing correlation between immune genes and feature genes. (Figure 4l,m,n)

21.GSVA-NS.R: GSVA for train set. (Figure 6a)

22.GSVA-NS-test.R :GSVA for test set. (Figure 6b)

23.feature genes_NIR genes_cor.R: Correlation analysis between feature genes and neuroinflammatory response genes in train set. (Figure 6c)

24.feature genes_testNIR genes_cor.R: Correlation analysis between feature genes and neuroinflammatory response genes in test set. (Figure 6d)

25.cibersort_NIRCor.R: Correlation analysis between immune cells and neuroinflammatory response genes in train and test set. (Figure 6e,f)

26. enet.R: Construction and validation (internal and external) of Enet model. (part of Figure 8a,b)

26. KNN.R: Construction and validation (internal and external) of KNN model. (part of Figure 8a,b)

26. lightgbm.R: Construction and validation (internal and external) of LighGBM model. (part of Figure 8a,b)

26. logistic.R: Construction and validation (internal and external) of LR model. (part of Figure 8a,b)

26.XGboost-SHAP.R: Construction and validation (internal and external) of XGBoost model (part of Figure 8a,b), SHAP explanation of XGBoost model (Figure 9a,b,c,d,e,f,g), and five-fold validation for XGBoost model (Figure 8g).

27.model_comparism.R: Comparison between models through ROC curves and DCA curves (Figure 8c,d,e,f) and learning curves of all models (Figure 8g).

28.XGboost-SHAP_contain immune cells.R: SHAP interaction plot containing immune cells and feature genes (Figure 9h).

29.app.R: Online predicting website based on XGBoost mode.
