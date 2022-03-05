# Metagene Projection Methodology R script to run the analysis described in 
# Carmona-Saez et al. "Metagene projection characterizes GEN2.2 and CAL-1 as relevant human plasmacytoid den-dritic cell models"

# This code was adapted from R scripts developed by Pablo Tamayo and provided from the paper:
# Metagene projection for cross-platform, cross-species characterization of global transcriptional states
# P. Tamayo, D. Scanfeld, B. L. Ebert, M. A. Gillette, C. W. M. Roberts, and J.P. Mesirov  
# Proc. Natl. Acad. Sci. USA, 104: 5959-5964 2007. http://www.pnas.org/cgi/content/abstract/0701068104v1 

#
# To run: change the path and cut and paste (or source) this code inside the R GUI console
#
# While running this script calls "MetaGene.Projection" (see below) which will produce the following output files
# under the directory specified by "output.dir":
#
# Main output files:
# <identifier>.<date>_<time>.log.txt = File containing the parameters used in the run and the data ans time
# <identifier>.model_dataset.H.gct = projected model dataset
# <identifier>.all.H.cls =  projection of model + test datasets (cls phenotypes)
# <identifier>.all.H.gct = projection of model + test datasets (gct dataset)
# <identifier>.heatmap.pdf = heat map of projection
# <identifier>.heatmap.sorted.pdf = heat map of projection sorted inside each phenotype
# <identifier>.2D.proj.pdf = 2D biplot of projected model and test datasets
# <identifier>.H.htree.pdf = hierarchical tree built on the projected model and test datasets
# <identifier>.pred.gct = projection-based SVM prediction results (gct dataset)
# <identifier>.pred.txt = projection-based SVM prediction results (text file)
# <identifier>.H.mem.txt = clustering membership based on metagene with largest amplitude
#
# Other complementary output files:
# <identifier>.model.H.gct = H matrix from the NMF decomposition
# <identifier>.model.W.gct = W matrix from the NMF decomposition
# <identifier>.model_set.2.cls = model dataset after pre-preprocessing and refinement (cls phenotypes)
# <identifier>.model_set.2.gct = model dataset after pre-preprocessing and refinement (gct file)
# <identifier>.model_set.1.cls = model dataset after pre-preprocessing and before refinement (cls phenotypes)
# <identifier>.model_set.1.gct = model dataset after pre-preprocessing and before refinement (gct files)
# <identifier>.model_set.0.cls = model dataset before pre-preprocessing but containing samples after refinement (cls phenotypes)
# <identifier>.model_set.0.gct = model dataset before pre-preprocessing but containing samples after refinement (gct file)
# <identifier>.htree.pdf = hierarchical tree on original pre-projection dataset
# <identifier>.all.cls = consolidated model + test dataset in the space of common genes (cls phenotypes)
# <identifier>.all.gct = consolidated model + test dataset in the space of common genes (gct dataset)
# <identifier>.prelim.pred.txt = preliminary projection-based SVM prediction results (used in refinement) (text file)

MP.library.location  <-  "A:/CarmonaSaezMetagenes/additionalFile5/MP.Library.R"
source(MP.library.location, verbose=T, max.deparse.length=9999)   # Load Metagene Projection library 

# Define model & test datasets and parameters

model.dataset.table <- # Defines the input model dataset and pre-processing options for it
  list( # GSE28490
     gct.file = "A:/CarmonaSaezMetagenes/additionalFile5/GSE28490.gct", 
     cls.file = "A:/CarmonaSaezMetagenes/additionalFile5/GSE28490.cls",
     column.subset="ALL",
     column.sel.type = "samples",  # Selection type: "sample": or "phenotypes"
     thres = "NULL",               # Threshold to apply to dataset before projection
     ceil = "NULL",                # Ceiling to apply to dataset before projection
     fold = "NULL",                     # Fold change (max/min) for variation filter before projection
     delta = "NULL",                  # Absolute difference (max - min) for variation filter before projection
     norm = 6)                     # Normalization before projection (default 6 column-rank and rescaling normalization)

test.datasets.table <-  # Defines one or more input test datasets and pre-processing optio  ns for each of them
   list(
      list( # Gen2.2
         gct.file = "A:/CarmonaSaezMetagenes/additionalFile5/Gen2.2.gct", 
         cls.file = "A:/CarmonaSaezMetagenes/additionalFile5/Gen2.2.cls", 
         column.subset = "ALL",
         column.sel.type = "samples",  
         thres = "NULL",
         ceil = "NULL",
         norm = 6), 
      list(  # GSE28491 
         gct.file = "A:/CarmonaSaezMetagenes/additionalFile5/GSE28491.gct", 
         cls.file = "A:/CarmonaSaezMetagenes/additionalFile5/GSE28491.cls",
         column.subset =  "ALL", 
         column.sel.type = "samples",  
         thres = "NULL",
         ceil = "NULL",
         norm = 6),
      list(  # GSE12507
         gct.file = "A:/CarmonaSaezMetagenes/additionalFile5/GSE12507.gct", 
         cls.file = "A:/CarmonaSaezMetagenes/additionalFile5/GSE12507.cls",
         column.subset ="ALL", 
         column.sel.type = "samples",  
         thres = "NULL",
         ceil = "NULL",
         norm = 6),
      list(  # GSE35457
         gct.file = "A:/CarmonaSaezMetagenes/additionalFile5/GSE35457.gct", 
         cls.file = "A:/CarmonaSaezMetagenes/additionalFile5/GSE35457.cls",
         column.subset = "ALL",
         column.sel.type = "samples",  
         thres = "NULL",
         ceil = "NULL",
         norm = 6)
     )

# Define parameters for this specific run (see detailed definitions below)

identifier          <-    "Gen2.2"
k.proj              <-           9
alg                 <-    "NMF.div"
niter               <-        2000
seed                <-        NULL           #Random seed to initialize metagene matrices. It can be set to a given value to get same results across factorizations
nchar.phen          <-           2
postprojnorm        <-        TRUE
use.biplot          <-        TRUE
heatmap.row.norm    <-       FALSE
heatmap.cmap.type   <-           6
use.feature.names   <-       FALSE
high.conf.thres     <-         0.3
output.dir          <-     "A:/CarmonaSaezMetagenes/additionalFile5/Results/"
kernel              <-     "radial"
cost                <-           1
gamma               <-           5
theta               <-        0.05
model.set.refinement <-          T

# These are the symbols and colors to use for each phenotype in the model and test sets 
#          model samples:   square symbols
#                  color         symbol      phenotype
legend.list <- c("blue",          22,        # B-cells
                 "red",           22,        # CD4-T-cells
                 "green",         22,        # CD8-T-cells
                 "steelblue",     22,        # eosinophils
                 "chocolate",     22,        # mDCs
                 "orchid",        22,        # monocytes
                 "yellow",        22,        # neutrophils
                 "paleturquoise", 22,        # NK-cells
                 "darkred",       22,        # pDCs
        #       test samples:    circles and other symbols
                 "black",         8,        # Gen2.2
                 "orchid",        23,        # CD14-monocyte
                 "blue",          23,        # CD19-B-cell
                 "red",           23,        # CD4 T cell
                 "green",         23,        # CD8 T cell
                 "steelblue",     23,        # eosinophil
                 "yellow",        23,        # neutrophil
                 "paleturquoise", 23,        # NK-cell
                 "gray60",        24,        # CAL1
                 "darkred",       25,        # pDC
                 "greenyellow",   25,        # CD14-CD16
                 "darksalmon",    25,        # CD14
                 "gold2",         25,        # CD141
                 "lightblue1",    25,        # CD16
                 "hotpink",       25         # CD1c
)

symbol.scaling <- 0.55
col <- legend.list[seq(1, length(legend.list), 2)]
symbs <- as.numeric(legend.list[seq(2, length(legend.list), 2)])

# This is the call to the Metagene Projection function:

MetaGene.Projection(                           # Runs the entire methodology
                                               #   (except for the GSEA analysis of metagenes and the model (k) selection)
  model.dataset.table = model.dataset.table,   # R list with model dataset parameters (see model.dataset.table above)
  test.datasets.table = test.datasets.table,   # R list with test dataset(s) parameters (see model.dataset.table above)
  identifier = identifier,                     # Prefix to name output files
  k.proj = k.proj,                             # Number of metagenes in projection
  alg = alg,                                   # Algorithm for Metagene Projection (default NMF.div):
                                               #    "NMF.div" : Non-Negative Matrix Factorization using the divergence cost
                                               #  (other algorithms for projection are internally supported but have not being tested)
  niter = niter,                               # Number of algorithm iterations (default: 2000)
  seed = seed,                                 # Random seed to initialize metagene matrices (default: 1234)
  nchar.phen =  nchar.phen,                    # Number of characters to use to identify classes from the CLS files
  postprojnorm = postprojnorm,                 # TRUE or FALSE: apply post-projection normalization (i.e. scale points to unit
                                               #     hypersphere, default: T)
  heatmap.row.norm = heatmap.row.norm,         # TRUE or FALSE: row-normalize (standardize) the rows in the heat map (default F)
  heatmap.cmap.type = heatmap.cmap.type,       # 1 = vintage pinkogram, 2 = scale of grays, 4 = high-resolution pinkogram,
                                               #   6 = redish color map for metagene factors (default: 6)
  high.conf.thres = high.conf.thres,           # Confidence threshold (Brier score) to separate call from no-calls (default 0.3)
  output.dir = output.dir,                     # Output directory where the resureflting output files will be produced
  col = col,                                   # Colors for the legend symbols for phenotypes: first model and then test dataset(s)
  symbs = symbs,                               # Plotting symbols for phenotypes: first model and then test dataset(s)
  symbol.scaling = symbol.scaling,             # Graphical scaling for symbols in plots and plot legends (default: 1)
  kernel = kernel,                             # Kernel function for SVM: "radial" or "linear" (default: "radial")
  cost = cost,                                 # Cost parameter for SVM (default: 1)
  gamma = gamma,                               # Gamma coefficient for radial base function kernel (default:  0.05 )
  model.set.refinement = model.set.refinement) # TRUE or FALSE: perform model set refinement (default: T)

# end of metagene projection script example

