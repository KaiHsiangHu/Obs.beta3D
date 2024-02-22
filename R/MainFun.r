#' Observed beta diversity with order q for TD, PD, and FD
#' 
#' \code{Obsbeta3D}: compute observed 3D alpha, beta, gamma diversity as well as dissimilarity indices.
#' 
#' @param data (a) For \code{datatype = "abundance"}, species abundance data for a single dataset can be input as a \code{matrix/data.frame} (species-by-assemblage); data for multiple datasets can be input as a \code{list} of \code{matrices/data.frames}, with each matrix representing a species-by-assemblage abundance matrix for one of the datasets.\cr
#' (b) For \code{datatype = "incidence_raw"}, data for a single dataset with N assemblages can be input as a \code{list} of \code{matrices/data.frames}, with each matrix representing a species-by-sampling-unit incidence matrix for one of the assemblages; data for multiple datasets can be input as multiple lists.
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' @param q a numerical vector specifying the diversity orders. Default is \code{c(0, 1, 2)}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence/occurrence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Set \code{nboot = 0} to skip the bootstrap procedures. Default is \code{nboot = 10}. If more accurate results are required, set \code{nboot = 100} (or \code{nboot = 200}).
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param PDtree (required argument for \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled assemblage. 
#' @param PDreftime (argument only for \code{diversity = "PD"}), a numerical value specifying reference time for PD. Default is \code{PDreftime = NULL} (i.e., the age of the root of \code{PDtree}).  
#' @param PDtype (argument only for \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"} (effective total branch length) or \code{PDtype = "meanPD"} (effective number of equally divergent lineages). Default is \code{PDtype = "meanPD"}, where \code{meanPD = PD/tree depth}.
#' @param FDdistM (required argument for \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled dataset. 
#' @param FDtype (argument only for \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_value"} for FD under a specified threshold value, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{FDtype = "AUC"}.  
#' @param FDtau (argument only for \code{diversity = "FD"} and \code{FDtype = "tau_value"}), a numerical value between 0 and
#'  1 specifying the tau value (threshold level) that will be used to compute FD. If \code{FDtype = NULL} (default), 
#'  then threshold level is set to be the mean distance between any two individuals randomly selected from the pooled 
#'  dataset (i.e., quadratic entropy). 
#' @param FDcut_number (argument only for \code{diversity = "FD"} and \code{FDtype = "AUC"}), a numeric number to cut [0, 1] interval into equal-spaced sub-intervals to obtain the AUC value by integrating the tau-profile. Equivalently, the number of tau values that will be considered to compute the integrated AUC value. Default is \code{FDcut_number = 30}. A larger value can be set to obtain more accurate AUC value.
#' 
#' @import magrittr
#' @import ggplot2
#' @import abind
#' @import iNEXT.3D
#' @import future.apply
#' @import tibble
#' @import dplyr
#' @import tidytree
#' @importFrom tidyr gather
#' @importFrom phyclust get.rooted.tree.height
#' @importFrom stats rmultinom
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom stats optimize
#' 
#' @return Return a list of several data frames including three diversity (gamma, alpha, and beta
#'  diversity) and four dissimilarity measures.\cr 
#'  
#'  The output in each data frame includes: 
#'  \item{Dataset}{the name of dataset.}
#'  \item{Order.q}{the diversity order of q.} 
#'  \item{Size}{the sample size.}
#'  \item{Value}{the observed gamma, alpha, beta diversity and dissimilarity value.}
#'  \item{s.e.}{standard error of standardized estimate.}
#'  \item{LCL, UCL}{the bootstrap lower and upper confidence limits for the diversity/dissimilarity with a default significance level of 0.95.}
#'  \item{Diversity}{'TD' = 'Taxonomic diversity', 'PD' = 'Phylogenetic diversity', 'meanPD' = 'Mean phylogenetic diversity', 'FD_tau' = 'Functional diversity (given tau)', 'FD_AUC' = 'Functional diversity (AUC)'}
#'  \item{Type}{gamma, alpha, beta diversity or dissimilarity measures.}  
#'  \item{Reftime}{the reference time for PD.}
#'  \item{Tau}{the threshold of functional distinctiveness between any two species for FD (under \code{FDtype = "tau_value"}).}
#'  
#'  
#' 
#' @examples
#' ## Taxonomic diversity for abundance data
#' data(beetle_abu)
#' output1 = Obsbeta3D(data = beetle_abu, diversity = 'TD', datatype = 'abundance', 
#'                     nboot = 10, conf = 0.95)
#' output1
#' 
#' 
#' ## Phylogenetic diversity for abundance data
#' data(beetle_abu)
#' data(beetle_tree)
#' output2 = Obsbeta3D(data = beetle_abu, diversity = 'PD', datatype = 'abundance', 
#'                     nboot = 10, conf = 0.95,
#'                     PDtree = beetle_tree, PDreftime = NULL)
#' output2
#' 
#' 
#' ## Functional diversity for abundance data under single threshold
#' data(beetle_abu)
#' data(beetle_distM)
#' output3 = Obsbeta3D(data = beetle_abu, diversity = 'FD', datatype = 'abundance', 
#'                     nboot = 10, conf = 0.95, 
#'                     FDdistM = beetle_distM, FDtype = 'tau_value', FDtau = NULL)
#' output3
#' 
#' 
#' ## Functional diversity for abundance data with thresholds integrating from 0 to 1
#' data(beetle_abu)
#' data(beetle_distM)
#' output4 = Obsbeta3D(data = beetle_abu, diversity = 'FD', datatype = 'abundance', 
#'                     nboot = 10, conf = 0.95, 
#'                     FDdistM = beetle_distM, FDtype = 'AUC', FDcut_number = 30)
#' output4
#' 
#' 
#' ## Taxonomic diversity for incidence data
#' data(beetle_inc)
#' output5 = Obsbeta3D(data = beetle_inc, diversity = 'TD', datatype = 'incidence_raw', 
#'                     nboot = 10, conf = 0.95)
#' output5
#' 
#' 
#' ## Phylogenetic diversity for incidence data
#' data(beetle_inc)
#' data(beetle_tree)
#' output6 = Obsbeta3D(data = beetle_inc, diversity = 'PD', datatype = 'incidence_raw', 
#'                     nboot = 10, conf = 0.95, 
#'                     PDtree = beetle_tree, PDreftime = NULL, PDtype = 'PD')
#' output6
#' 
#' 
#' ## Functional diversity for incidence data under single threshold
#' data(beetle_inc)
#' data(beetle_distM)
#' output7 = Obsbeta3D(data = beetle_inc, diversity = 'FD', datatype = 'incidence_raw', 
#'                     nboot = 10, conf = 0.95, 
#'                     FDdistM = beetle_distM, FDtype = 'tau_value', FDtau = NULL)
#' output7
#' 
#' 
#' ## Functional diversity for incidence data with thresholds integrating from 0 to 1
#' data(beetle_inc)
#' data(beetle_distM)
#' output8 = Obsbeta3D(data = beetle_inc, diversity = 'FD', datatype = 'incidence_raw', 
#'                     nboot = 10, conf = 0.95, 
#'                     FDdistM = beetle_distM, FDtype = 'AUC', FDcut_number = 30)
#' output8
#' 
#' @export
Obsbeta3D = function(data, diversity = 'TD', q = seq(0, 2, 0.25), datatype = 'abundance',
                     nboot = 20, conf = 0.95, PDtree = NULL, PDreftime = NULL, PDtype = 'meanPD',
                     FDdistM = NULL, FDtype = 'AUC', FDtau = NULL, FDcut_number = 50) {
  
  ## Check parameter setting
  if (is.na(pmatch(diversity, c("TD", "PD", "FD")))) stop("invalid diversity")
  
  if (!inherits(q, "numeric"))
    stop("invlid class of order q, q should be a postive value/vector of numeric object", call. = FALSE)
  if (min(q) < 0){
    warning("ambigous of order q, we only compute postive q.", call. = FALSE)
    q <- q[q >= 0]
  }
  
  if (datatype == "incidence" | datatype == "incidence_freq") stop('Please try datatype = "incidence_raw".')  
  if(is.na(pmatch(datatype, c("abundance", "incidence_raw"))))
    stop("invalid datatype")
  
  if ((nboot < 0) | (is.numeric(nboot) == F)) stop('Please enter non-negative integer for nboot.', call. = FALSE)
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf) == F)) stop('Please enter value between zero and one for confident interval.', call. = FALSE)
  
  if (! (is.null(PDreftime) | inherits(PDreftime, "numeric")))
    stop("invalid class of reference time, PDreftime should be a postive value of numeric object.", call. = FALSE)
  if (length(PDreftime) > 1)
    stop("PDreftime can only accept a value instead of a vector.", call. = FALSE)
  
  if(is.na(pmatch(PDtype, c("PD", "meanPD"))))
    stop("Incorrect type of phylogenetic diversity type, please use either 'PD' or 'meanPD'.", call. = FALSE)
  
  if (FDtype == "tau_values") stop('Please try FDtype = "tau_value".')  
  if(is.na(pmatch(FDtype, c("AUC", "tau_value"))))
    stop("Incorrect type of functional diversity type, please use either 'AUC' or 'tau_value'", call. = FALSE)
  
  if (! (is.null(FDtau) | inherits(FDtau, "numeric")))
    stop("invalid class of tau value, FDtau should be a postive value between zero and one.", call. = FALSE)
  if (length(FDtau) > 1)
    stop("FDtau only accept a value instead of a vector.", call. = FALSE)
  
  if (!inherits(FDcut_number, "numeric"))
    stop("invalid class of FD cut number, FDcut_number should be a postive value.", call. = FALSE)
  if (FDcut_number < 2)
    stop("invalid FDcut_number, FDcut_number should be a postive value larger than one.", call. = FALSE)
  if (length(FDcut_number) > 1)
    stop("FDcut_number only accept a value instead of a vector.", call. = FALSE)
  
  ##
  
  if (datatype == 'abundance') {
    
    if ( inherits(data, "data.frame") | inherits(data, "matrix") ) data = list(Dataset_1 = data)
    
    if ( inherits(data, "list") ){
      
      if (is.null(names(data))) dataset_names = paste0("Dataset_", 1:length(data)) else dataset_names = names(data)
      Ns = sapply(data, ncol)
      data_list = data
      
    }
    
  }
  
  if (datatype == 'incidence_raw') {
    
    if (!inherits(data, "list"))
      stop("Invalid data format for incidence raw data. Please refer to example of iNEXTbeta3D.", call. = FALSE)
    
    if ( inherits(data, "list") & (inherits(data[[1]], "data.frame") | inherits(data[[1]], "matrix")) ) data = list(Dataset_1 = data)
    
    if (! ( inherits(data[[1]][[1]], "data.frame") | inherits(data[[1]][[1]], "matrix") ) )
      stop("Invalid data format for incidence raw data. Please refer to example of iNEXTbeta3D.", call. = FALSE)
    
    if ( sum( sapply(1:length(data), function(i) ( length(unique(sapply(data[[i]], nrow))) != 1 | 
                                                   length(unique(sapply(data[[i]], ncol))) != 1 ) ) ) > 0 )
      stop("Number of species (row) or sampling units (column) should be the same within each dataset. Please check you data or refer to example of iNEXTbeta3D.", call. = FALSE)
    
    data = lapply(data, function(x) {
      
      if (is.null(rownames(x[[1]]))) tmp = x else {
        
        nT = ncol(x[[1]])
        if (nT <= 3) stop("Number of sampling units of some datasets is too less. Please add more sampling units data.", call. = FALSE)
        
        tmp = lapply(x, function(i) data.frame(i) %>% rownames_to_column(var = "Species"))
        
        tmp1 = tmp[[1]]
        for (i in 2:length(tmp)) {
          tmp1 = full_join(tmp1, tmp[[i]], by = "Species")
        }
        tmp1 = tmp1 %>% column_to_rownames(var = "Species")
        
        if (nrow(tmp1) != nrow(x[[1]]))
          stop("Species names (rownames) should be matched within each dataset. Please check you data or refer to example of iNEXTbeta3D.", call. = FALSE)
        
        if (sum(!as.matrix(tmp1) %in% c(0,1)) > 0)
          stop("The data for datatype = 'incidence_raw' can only contain values zero (undetected) or one (detected). Please transform values to zero or one.", call. = FALSE)
        
        tmp = lapply(1:length(tmp), function(i) tmp1[, ((i-1)*nT+1):(i*nT)])
        names(tmp) = if (is.null(data)) paste0("Assemblage_", 1:length(data)) else names(x)
      }
      
      return(tmp)
    })
    
    if (is.null(names(data))) dataset_names = paste0("Dataset_", 1:length(data)) else dataset_names = names(data)
    Ns = sapply(data, length)
    data_list = data
  }
  
  
  if (datatype == 'abundance') {
    
    pool.name <- lapply(data_list, function(x) rownames(x)) %>% unlist %>% unique
    
  } else if (datatype == 'incidence_raw') {
    
    pool.name <- lapply(data_list, function(x) lapply(x, function(y) rownames(y))) %>% unlist %>% unique
    
  }
  
  if (diversity == "PD") {
    
    if (sum(c(duplicated(PDtree$tip.label), duplicated(PDtree$node.label[PDtree$node.label!=""])))>0)
      stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
    
    if ( is.null(pool.name) )
      stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
    
    if (sum(pool.name %in% PDtree$tip.label) != length(pool.name))
      stop("Data and tree tip label contain unmatched species", call. = FALSE)
  }
  
  if (diversity == "FD") {
    
    if (is.null(rownames(FDdistM)))
      stop('The species names are not provided in distance matrix.', call. = FALSE)
    
    if( is.null(pool.name) )
      stop("Row names of data must be the species names that match row names in distance matrix and thus can not be empty.", call. = FALSE)
    
    if (sum(pool.name %in% rownames(FDdistM)) != length(pool.name))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  ##
  
  
  if (is.null(conf)) conf = 0.95
  tmp = qnorm(1 - (1 - conf)/2)
  
  if (diversity == 'FD' & FDtype == 'tau_value' & is.null(FDtau) == T) {
    
    if (datatype == 'abundance') {
      
      tt <- lapply(data_list, rowSums)
      tt = lapply(tt, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
      pdata = tt[[1]]
      
      if (length(tt) > 1) {
        for(i in 2:length(tt)){
          pdata = full_join(pdata, tt[[i]], by = "Species")
        }
      }
      
      pdata[is.na(pdata)] = 0
      pdata = pdata %>% column_to_rownames("Species")
      pdata = rowSums(pdata)
      
      order_sp <- match(names(pdata),rownames(FDdistM))
      FDdistM <- FDdistM[order_sp,order_sp]
      pdata <- matrix(pdata/sum(pdata), ncol = 1)
      
    } else if (datatype == 'incidence_raw') {
      
      tt <- lapply(data_list, function(x) {tt = Reduce('+', x); tt[tt > 1] = 1; rowSums(tt) })
      tt = lapply(tt, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
      pdata = tt[[1]]
      
      if (length(tt) > 1) {
        for(i in 2:length(tt)){
          pdata = full_join(pdata, tt[[i]], by = "Species")
        }
      }
      
      pdata[is.na(pdata)] = 0
      pdata = pdata %>% column_to_rownames("Species")
      pdata = rowSums(pdata)
      
      order_sp <- match(names(pdata),rownames(FDdistM))
      FDdistM <- FDdistM[order_sp,order_sp]
      pdata <- matrix(pdata/sum(pdata), ncol = 1)
      
    }
    
    FDtau <- sum ( (pdata %*% t(pdata) ) * FDdistM) # dmean
  }
  
  if (diversity == 'PD') {
    
    if (datatype == "abundance") 
      
      if (length(data_list) > 1) {
        
        pool.data = data_list[[1]] %>% data.frame %>% rownames_to_column()
        for (i in 2:length(data_list)) 
          pool.data = full_join(pool.data, data_list[[i]] %>% data.frame %>% rownames_to_column(), 'rowname')
        pool.data[is.na(pool.data)] = 0
        pool.data = pool.data %>% column_to_rownames() %>% rowSums
        
      } else pool.data = do.call(cbind, data_list) %>% rowSums
      
      if (datatype == 'incidence_raw') pool.data = do.call(cbind,lapply(data_list, function(x) do.call(cbind,x)) ) %>% rowSums
      
      pool.name = names(pool.data[pool.data>0])
      tip = PDtree$tip.label[-match(pool.name, PDtree$tip.label)]
      mytree = ape::drop.tip(PDtree, tip)
      
      # H_max = get.rooted.tree.height(mytree)
      H_max = max(ape::node.depth.edgelength(mytree))
      
      if(is.null(PDreftime)) { reft = H_max
      } else if (PDreftime <= 0) { stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.", call. = FALSE)
      } else { reft = PDreftime }
      
  }
  
  for_each_dataset = function(data, dataset_name, N) {
    
    #data
    if (datatype == 'abundance') {
      
      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector
      
    }
    
    if (datatype == 'incidence_raw') {
      
      sampling_units = sapply(data, ncol)
      if (length(unique(sampling_units)) > 1) stop("unsupported data structure: the sampling units of all datasets must be the same.")
      if (length(unique(sampling_units)) == 1) n = unique(sampling_units)
      
      gamma = Reduce('+', data)
      gamma[gamma>1] = 1
      data_gamma_raw = gamma
      data_gamma_freq = c(n, rowSums(gamma))
      
      data_alpha_freq = sapply(data, rowSums) %>% c(n, .)
      
      data_2D = apply(sapply(data, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
      
    }
    
    
    
    if (diversity == 'TD') {
      
      if (datatype == 'abundance') {
        
        gamma = ObsAsy3D(as.numeric(data_gamma), 'TD', q = q, datatype = "abundance", nboot = 0, method = 'Observed')
        
        alpha = ObsAsy3D(as.numeric(data_alpha), 'TD', q = q, datatype = "abundance", nboot = 0, method = 'Observed')
        
      }
      
      if (datatype == 'incidence_raw') {
        
        gamma = ObsAsy3D(as.numeric(data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0, method = 'Observed')
        
        alpha = ObsAsy3D(as.numeric(data_alpha_freq), diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0, method = 'Observed')
        
        
      }
      
      gamma = gamma[,c(3,2)] %>% set_colnames(c('Value', 'Order.q'))
      
      alpha = alpha[,c(3,2)] %>% set_colnames(c('Value', 'Order.q'))
      
      alpha$Value = alpha$Value / N
      
      beta = alpha
      beta$Value = gamma$Value/alpha$Value
      C = beta %>% mutate(Value = ifelse(Order.q==1, log(Value)/log(N), (Value^(1-Order.q) - 1)/(N^(1-Order.q)-1)))
      U = beta %>% mutate(Value = ifelse(Order.q==1, log(Value)/log(N), (Value^(Order.q-1) - 1)/(N^(Order.q-1)-1)))
      V = beta %>% mutate(Value = (Value-1)/(N-1))
      S = beta %>% mutate(Value = (1/Value-1)/(1/N-1))
      
      if(nboot>1){
        
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = bootstrap_population[,k]))
            
            bootstrap_data_gamma = rowSums(bootstrap_sample)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma > 0]
            bootstrap_data_alpha = as.matrix(bootstrap_sample) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha > 0]
            
            gamma = ObsAsy3D(as.numeric(bootstrap_data_gamma), diversity = 'TD', q = q, datatype = "abundance", nboot = 0, method = 'Observed')
            
            alpha = ObsAsy3D(as.numeric(bootstrap_data_alpha), diversity = 'TD', q = q, datatype = "abundance", nboot = 0, method = 'Observed')
            
          }
          
          if (datatype == 'incidence_raw') {
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            
            raw = lapply(1:ncol(bootstrap_population), function(j){
              
              lapply(1:nrow(bootstrap_population), function(i) rbinom(n = n, size = 1, prob = bootstrap_population[i,j])) %>% do.call(rbind,.)
              
            })
            
            gamma = Reduce('+', raw)
            gamma[gamma > 1] = 1
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))
            
            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)
            
            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]
            
            gamma = ObsAsy3D(bootstrap_data_gamma_freq, diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0, method = 'Observed')
            
            alpha = ObsAsy3D(bootstrap_data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0, method = 'Observed')
            
          }
          
          gamma = gamma$qTD
          
          alpha = alpha$qTD
          alpha = alpha / N
          
          beta = data.frame(Value = gamma/alpha, q)
          
          C = (beta %>% mutate(Value = ifelse(q == 1,log(Value)/log(N),(Value^(1 - q) - 1)/(N^(1 - q) - 1))))$Value
          U = (beta %>% mutate(Value = ifelse(q == 1,log(Value)/log(N),(Value^(q - 1) - 1)/(N^(q - 1) - 1))))$Value
          V = (beta %>% mutate(Value = (Value - 1)/(N - 1)))$Value
          S = (beta %>% mutate(Value = (1/Value - 1)/(1/N - 1)))$Value
          
          beta = beta$Value
          
          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          
        }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
        
      } else {
        
        se = matrix(0, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)
        
      }
      
    }
    
    if (diversity == 'PD') {
      
      if (datatype == 'abundance') {
        
        aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = reft)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        
        gamma = iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), q = q, reft = reft, cal = "PD", nt = n) %>% 
          matrix(., ncol = 1) %>% data.frame() %>% cbind(q) %>% set_colnames(c('Value', 'Order.q'))
        
        
        aL_table_alpha = c()
        
        for (i in 1:N){
          
          x = data[data[,i]>0,i]
          names(x) = rownames(data)[data[,i]>0]
          
          aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          aL_table_alpha = rbind(aL_table_alpha, aL_table)
          
        }
        
        
        qPDm = iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, reft = reft, cal = "PD", nt = n) %>%
          matrix(., ncol = 1)
        qPDm = qPDm/N
        alpha = qPDm %>% data.frame() %>% cbind(q) %>% set_colnames(c('Value', 'Order.q'))
        
      }
      
      if (datatype == 'incidence_raw') {
        
        aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = as.matrix(data_gamma_raw), datatype = "incidence_raw", refT = reft, rootExtend = T)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        gamma = iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis=as.matrix(aL_table_gamma$branch.length), q = q, reft = reft, cal = "PD", nt = n) %>% 
          matrix(., ncol = 1) %>% data.frame() %>% cbind(q) %>% set_colnames(c('Value', 'Order.q'))
        
        aL_table_alpha = c()
        
        for (i in 1:N){
          
          x = data[[i]]
          
          aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          aL_table_alpha = rbind(aL_table_alpha, aL_table)
          
        }
        
        alpha = (iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, reft = reft, cal = "PD", nt = n)/N) %>% 
          matrix(., ncol = 1) %>% cbind(q) %>% data.frame() %>% set_colnames(c('Value', 'Order.q'))
        
        
      }
      
      if (PDtype == 'meanPD') {
        gamma$Value = gamma$Value/reft
        alpha$Value = alpha$Value/reft
      }
      
      beta = alpha
      beta$Value = gamma$Value/alpha$Value
      
      C = beta %>% mutate(Value = ifelse(Order.q == 1, log(Value)/log(N), (Value^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1)))
      U = beta %>% mutate(Value = ifelse(Order.q == 1, log(Value)/log(N), (Value^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1)))
      V = beta %>% mutate(Value = (Value - 1)/(N - 1))
      S = beta %>% mutate(Value = (1/Value - 1)/(1/N - 1))
      
      if(nboot>1){
        
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            tree_bt = PDtree
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data)),] %>% matrix(ncol = ncol(data))
            
            if ( nrow(p_bt) > nrow(data) & sum(unseen_p) > 0 ){
              
              unseen = unseen_p[which(rowSums(unseen_p) > 0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data)),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data), unseen_name)
              
              bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              x_bt = bootstrap_sample
              
              rownames(x_bt) = rownames(p_bt)
              
              if ( sum(x_bt[-(1:nrow(data)),])>0 ){
                
                g0_hat = apply(data, 2, function(x){
                  
                  n = sum(x)
                  f1 = sum(x == 1)
                  f2 = sum(x == 2)
                  
                  aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
                  
                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  g1 = aL$branch.length[aL$branch.abun == 1] %>% sum
                  g2 = aL$branch.length[aL$branch.abun == 2] %>% sum
                  g0_hat = ifelse( g2 > ((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
                  if(is.na(g0_hat)) {g0_hat <- 0 }
                  g0_hat
                  
                })
                
                te = (x_bt[1:nrow(data),]*(data == 0))>0
                used_length = sapply(1:ncol(data), function(i) { 
                  
                  if (sum(te[,i]) == 0) return(0) else {
                    
                    iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x_bt[1:nrow(data),i], rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum
                    
                  }
                  
                })
                
                g0_hat = g0_hat - used_length
                g0_hat[g0_hat < 0] = 0
                
                unseen_sample = x_bt[-(1:nrow(data)),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = ncol(x_bt), byrow = T)
                
                L0_hat = sapply(1:length(g0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (g0_hat[i] / nrow(unseen)) else 0 )
                
                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 0)] = 0
                
                for (i in 1:length(L0_hat)){
                  
                  tip = list(edge = matrix(c(2,1),1,2),
                             tip.label = unseen_name[i],
                             edge.length = L0_hat[i],
                             Nnode = 1)
                  class(tip) = "phylo"
                  
                  tree_bt = tree_bt + tip
                  
                }
                
              } else {
                
                x_bt = x_bt[1:nrow(data),]
                p_bt = p_bt[1:nrow(data),]
                
              }
              
            } else {
              
              p_bt = p_bt[1:nrow(data),]
              x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              rownames(x_bt) = rownames(data)
              
            }
            
            bootstrap_data_gamma = rowSums(x_bt)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]
            bootstrap_data_alpha = as.matrix(x_bt) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]
            
            aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = bootstrap_data_gamma, rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            gamma = as.vector(iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), q = q, reft = reft, cal = "PD", nt = n))
            
            
            aL_table_alpha = c()
            
            for (i in 1:N){
              
              x = x_bt[,i]
              names(x) = rownames(p_bt)
              x = x[x_bt[,i]>0]
              
              aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = x, rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
              
              aL_table_alpha = rbind(aL_table_alpha, aL_table)
              
            }
            
            alpha = as.vector((iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, reft = reft, cal = "PD", nt = n)/N))
            
          }
          
          if (datatype == 'incidence_raw') {
            
            tree_bt = PDtree
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data[[1]])),] %>% matrix(ncol=N)
            
            if ( nrow(p_bt) > nrow(data[[1]]) & sum(unseen_p)>0 ){
              
              unseen = unseen_p[which(rowSums(unseen_p)>0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data[[1]])),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data[[1]]), unseen_name)
              
              raw = lapply(1:ncol(p_bt), function(j){
                
                lapply(1:nrow(p_bt), function(i) rbinom(n=n, size=1, prob=p_bt[i,j])) %>% do.call(rbind,.)
                
              })
              
              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
              
              if ( lapply(1:length(raw), function(i) raw[[i]][-(1:nrow(data[[1]])),]) %>% do.call(sum,.)>0 ){
                
                R0_hat = sapply(data, function(x){
                  
                  nT = ncol(x)
                  Q1 = sum(rowSums(x)==1)
                  Q2 = sum(rowSums(x)==2)
                  
                  aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
                  
                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  R1 = aL$branch.length[aL$branch.abun == 1] %>% sum
                  R2 = aL$branch.length[aL$branch.abun == 2] %>% sum
                  R0_hat = ifelse( R2>((R1*Q2)/(2*Q1)) , ((nT-1)/nT)*(R1^2/(2*R2)) , ((nT-1)/nT)*(R1*(Q1-1)/(2*(Q2+1))) )
                  if(is.na(R0_hat)) { R0_hat <- 0 }
                  R0_hat
                  
                })
                
                te = (sapply(raw, rowSums)[1:nrow(data[[1]]),]*(sapply(data, rowSums) == 0)) > 0
                used_length = sapply(1:N, function(i) {
                  
                  if (sum(te[,i]) == 0) return(0) else {
                    
                    iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = raw[[i]][1:nrow(data[[1]]),], datatype = "incidence_raw", rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum
                    
                  }
                  
                })
                
                R0_hat = R0_hat - used_length
                R0_hat[R0_hat < 0] = 0
                
                unseen_sample = sapply(raw, rowSums)[-(1:nrow(data[[1]])),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = N, byrow = T)
                
                L0_hat = sapply(1:length(R0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (R0_hat[i] / nrow(unseen)) else 0 )
                
                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 0)] = 0
                
                for (i in 1:length(L0_hat)){
                  
                  tip = list(edge = matrix(c(2,1), 1, 2),
                             tip.label = unseen_name[i],
                             edge.length = L0_hat[i],
                             Nnode = 1)
                  class(tip) = "phylo"
                  
                  tree_bt = tree_bt + tip
                  
                }
                
              } else raw = lapply(raw, function(i) i[1:nrow(data[[1]]),])
              
            } else {
              
              p_bt = p_bt[1:nrow(data[[1]]),]
              raw = lapply(1:ncol(p_bt), function(j){
                
                lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)
                
              })
              
              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
              
            }
            
            gamma = Reduce('+', raw)
            gamma[gamma > 1] = 1
            bootstrap_data_gamma_raw = gamma
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))
            
            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)
            
            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]
            
            aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = bootstrap_data_gamma_raw, datatype = "incidence_raw", rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            gamma = as.vector(iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), q = q, reft = reft, cal = "PD", nt = n))
            
            
            aL_table_alpha = c()
            
            for (i in 1:N){
              
              x = raw[[i]]
              
              aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
              
              aL_table_alpha = rbind(aL_table_alpha, aL_table)
              
            }
            
            alpha = as.vector((iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, reft = reft, cal = "PD", nt = n)/N))
            
          }
          
          if (PDtype == 'meanPD') {
            gamma = gamma/reft
            alpha = alpha/reft
          } 
          
          beta = data.frame(Value = gamma/alpha, Order.q = q)
          
          C = (beta %>% mutate(Value = ifelse(Order.q == 1, log(Value)/log(N), (Value^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Value
          U = (beta %>% mutate(Value = ifelse(Order.q == 1, log(Value)/log(N), (Value^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Value
          V = (beta %>% mutate(Value = (Value - 1)/(N - 1)))$Value
          S = (beta %>% mutate(Value = (1/Value - 1)/(1/N - 1)))$Value
          
          beta = beta$Value
          
          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          
        }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
        
      } else {
        
        se = matrix(0, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)
        
      }
      
    }
    
    if (diversity == 'FD') {
      
      FDdistM = as.matrix(FDdistM)
      
      FD_by_tau = function(data, distM, tau, datatype) {
        
        if (datatype == 'abundance') {
          
          zik = data
          zik = zik[rowSums(data)>0,]
          
          dij = distM
          dij = dij[rowSums(data)>0, rowSums(data)>0]
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            aik = (1 - dij/1) %*% as.matrix(zik)
          } else {
            dij[which(dij > tau, arr.ind = T)] = tau
            aik = (1 - dij/tau) %*% as.matrix(zik)
          }
          
          positive_id = rowSums(aik)>0
          gamma_x = rowSums(zik)[positive_id]
          gamma_a = rowSums(aik)[positive_id]
          gamma_v = gamma_x/gamma_a
          
          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
          
          gamma = iNEXT.3D:::FD_mle(ai_vi_gamma, q) %>% as.vector
          
          
          alpha_x = as.vector(as.matrix(zik))
          alpha_a = as.vector(aik)
          
          alpha_v = alpha_x/alpha_a
          alpha_v = rep(gamma_v,N)
          
          alpha_v = alpha_v[alpha_a>0]
          alpha_a = alpha_a[alpha_a>0]
          
          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
          
          alpha = (iNEXT.3D:::FD_mle(ai_vi_alpha, q)/N) %>% as.vector
          
        }
        
        if (datatype == 'incidence_raw') {
          
          data_gamma_freq = data$data_gamma_freq
          data_2D = data$data_2D
          
          gamma_Y = data_gamma_freq[-1]
          
          dij = distM
          dij = dij[gamma_Y > 0, gamma_Y > 0]
          gamma_Y = gamma_Y[gamma_Y > 0]
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            gamma_a = (1 - dij/1) %*% as.matrix(gamma_Y)
          } else {
            dij[which(dij > tau, arr.ind = T)] = tau
            gamma_a = (1 - dij/tau) %*% as.matrix(gamma_Y)
          }
          
          gamma_a[gamma_a > n] = n
          gamma_v = gamma_Y/gamma_a
          
          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
          
          gamma = iNEXT.3D:::FD_mle(ai_vi_gamma, q) %>% as.vector
          
          alpha_Y = data_2D[-1,]
          
          dij = distM
          dij = dij[rowSums(data_2D[-1,]) > 0, rowSums(data_2D[-1,])>0]
          alpha_Y = alpha_Y[rowSums(data_2D[-1,])>0,]
          
          dij[which(dij>tau, arr.ind = T)] = tau
          alpha_a = (1 - dij/tau) %*% as.matrix(alpha_Y)
          
          alpha_a[alpha_a > n] = n
          alpha_a = as.vector(alpha_a)
          
          alpha_v = rep(gamma_v, N)
          alpha_v = alpha_v[alpha_a > 0]
          alpha_a = alpha_a[alpha_a > 0]
          
          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
          
          alpha = (iNEXT.3D:::FD_mle(ai_vi_alpha, q)/N) %>% as.vector
          
        }
        
        return(data.frame(gamma,alpha))
        
      }
      
      if (FDtype == 'tau_value'){
        
        if (datatype == 'abundance') {
          
          output = FD_by_tau(data, FDdistM, FDtau, datatype = 'abundance')
          gamma = output$gamma
          alpha = output$alpha
          
          gamma = data.frame(gamma, q) %>% set_colnames(c('Value', 'Order.q'))
          
          alpha = data.frame(alpha, q) %>% set_colnames(c('Value', 'Order.q'))
          
          beta = alpha
          beta$Value = gamma$Value/alpha$Value
          
        }
        
        if (datatype == 'incidence_raw') {
          
          output = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, FDtau, datatype='incidence_raw')
          gamma = output$gamma
          alpha = output$alpha
          
          gamma = data.frame(gamma, q) %>% set_colnames(c('Value', 'Order.q'))
          
          alpha = data.frame(alpha, q) %>% set_colnames(c('Value', 'Order.q'))
          
          beta = alpha
          beta$Value = gamma$Value/alpha$Value
          
        }
        
      }
      
      if (FDtype == 'AUC'){
        
        cut = seq(0.00000001, 1, length.out = FDcut_number)
        width = diff(cut)
        
        if (datatype == 'abundance') {
          
          gamma_alpha_over_tau = lapply(cut, function(tau) {
            
            FD_by_tau(data, FDdistM, tau, datatype = 'abundance')
            
          })
          
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
          
          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
          
          gamma = colSums((left_limit + right_limit)/2)
          
          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
          
          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
          
          alpha = colSums((left_limit + right_limit)/2)
          
          beta_over_tau = gamma_over_tau/alpha_over_tau
          
          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
          
          beta = colSums((left_limit + right_limit)/2)
          
          gamma = data.frame(gamma, q) %>% set_colnames(c('Value', 'Order.q'))
          
          alpha = data.frame(alpha, q) %>% set_colnames(c('Value', 'Order.q'))
          
          beta = data.frame(beta, q) %>% set_colnames(c('Value', 'Order.q'))
          
        }
        
        if (datatype == 'incidence_raw') {
          
          gamma_alpha_over_tau = lapply(cut, function(tau) {
            
            FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, tau, datatype = 'incidence_raw')
            
          })
          
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
          
          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
          
          gamma = colSums((left_limit + right_limit)/2)
          
          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
          
          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
          
          alpha = colSums((left_limit + right_limit)/2)
          
          beta_over_tau = gamma_over_tau/alpha_over_tau
          
          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
          
          beta = colSums((left_limit + right_limit)/2)
          
          gamma = data.frame(gamma, q) %>% set_colnames(c('Value', 'Order.q'))
          
          alpha = data.frame(alpha, q) %>% set_colnames(c('Value', 'Order.q'))
          
          beta = data.frame(beta, q) %>% set_colnames(c('Value', 'Order.q'))
          
        }
        
      }
      
      C = beta %>% mutate(Value = ifelse(Order.q == 1, log(Value)/log(N), (Value^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1)))
      U = beta %>% mutate(Value = ifelse(Order.q == 1, log(Value)/log(N), (Value^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1)))
      V = beta %>% mutate(Value = (Value - 1)/(N - 1))
      S = beta %>% mutate(Value = (1/Value - 1)/(1/N - 1))
      
      if(nboot > 1){
        
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            f0_hat = nrow(p_bt) - nrow(data)
            
            distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), FDdistM, f0_hat, 'abundance')
            
            data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
            
            data_gamma = rowSums(data_bt)
            data_gamma = data_gamma[data_gamma>0]
            data_alpha = as.matrix(data_bt) %>% as.vector
            
            if (FDtype == 'tau_value'){
              
              output = FD_by_tau(data_bt, distance_matrix_bt, FDtau, datatype='abundance')
              gamma = output$gamma
              alpha = output$alpha
              beta=gamma/alpha
              
            }
            
            if (FDtype == 'AUC'){
              
              gamma_alpha_over_tau = lapply(cut, function(tau) {
                
                FD_by_tau(data_bt, distance_matrix_bt, tau, datatype = 'abundance')
                
              })
              
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              
              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
              
              gamma = colSums((left_limit + right_limit)/2)
              
              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              
              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
              
              alpha = colSums((left_limit + right_limit)/2)
              
              beta_over_tau = gamma_over_tau/alpha_over_tau
              
              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              
              beta = colSums((left_limit + right_limit)/2)
              
            }
            
          }
          
          if (datatype == 'incidence_raw') {
            
            p_bt = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            f0_hat = nrow(p_bt) - nrow(data_2D[-1,])
            
            distance_matrix_bt = Bootstrap_distance_matrix(c(n,rowSums(data_gamma_raw)), FDdistM, f0_hat, 'incidence_freq')
            
            raw = lapply(1:ncol(p_bt), function(j){
              
              lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)
              
            })
            
            gamma = Reduce('+', raw)
            gamma[gamma>1] = 1
            data_gamma_raw_bt = gamma
            data_gamma_freq_bt = c(n, rowSums(gamma))
            
            data_alpha_freq_bt = sapply(raw, rowSums) %>% c(n, .)
            
            data_2D_bt = apply(sapply(raw, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
            
            if (FDtype == 'tau_value'){
              
              output = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, FDtau, datatype = 'incidence_raw')
              gamma = output$gamma
              alpha = output$alpha
              beta = gamma/alpha
              
            }
            
            if (FDtype == 'AUC'){
              
              gamma_alpha_over_tau = lapply(cut, function(tau) {
                
                FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, datatype = 'incidence_raw')
                
              })
              
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              
              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
              
              gamma = colSums((left_limit + right_limit)/2)
              
              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              
              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
              
              alpha = colSums((left_limit + right_limit)/2)
              
              beta_over_tau = gamma_over_tau/alpha_over_tau
              
              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              
              beta = colSums((left_limit + right_limit)/2)
              
            }
            
          }
          
          beta = gamma/alpha
          
          beta = data.frame(Value = beta, Order.q = q)
          
          C = (beta %>% mutate(Value = ifelse(Order.q == 1, log(Value)/log(N), (Value^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Value
          U = (beta %>% mutate(Value = ifelse(Order.q == 1, log(Value)/log(N), (Value^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Value
          V = (beta %>% mutate(Value = (Value - 1)/(N - 1)))$Value
          S = (beta %>% mutate(Value = (1/Value - 1)/(1/N - 1)))$Value
          
          beta = beta$Value
          
          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          
          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
        
      } else {
        
        se = matrix(0, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)
        
      }
      
    }
    
    
    se = as.data.frame(se)
    
    if (diversity == "TD") index = "TD"
    if (diversity == "PD" & PDtype == "PD") index = "PD"
    if (diversity == "PD" & PDtype == "meanPD") index = "meanPD"
    if (diversity == "FD" & FDtype == "tau_value") index = "FD_tau"
    if (diversity == "FD" & FDtype == "AUC") index = "FD_AUC"
    
    gamma = gamma %>% mutate(Size = n,
                             s.e. = se$gamma,
                             LCL = Value - tmp * se$gamma,
                             UCL = Value + tmp * se$gamma,
                             Dataset = dataset_name,
                             Diversity = index) %>% select(c("Dataset", "Order.q", "Size", "Value", "s.e.", "LCL", "UCL", "Diversity"))
    
    alpha = alpha %>% mutate(Size = n,
                             s.e. = se$alpha,
                             LCL = Value - tmp * se$alpha,
                             UCL = Value + tmp * se$alpha,
                             Dataset = dataset_name,
                             Diversity = index) %>% select(c("Dataset", "Order.q", "Size", "Value", "s.e.", "LCL", "UCL", "Diversity"))
    
    beta = beta %>% mutate(  Size = n,
                             s.e. = se$beta,
                             LCL = Value - tmp * se$beta,
                             UCL = Value + tmp * se$beta,
                             Dataset = dataset_name,
                             Diversity = index) %>% select(c("Dataset", "Order.q", "Size", "Value", "s.e.", "LCL", "UCL", "Diversity"))
    
    C = C %>% mutate(        Size = n,
                             s.e. = se$C,
                             LCL = Value - tmp * se$C,
                             UCL = Value + tmp * se$C,
                             Dataset = dataset_name,
                             Diversity = index) %>% select(c("Dataset", "Order.q", "Size", "Value", "s.e.", "LCL", "UCL", "Diversity"))
    
    
    U = U %>% mutate(        Size = n,
                             s.e. = se$U,
                             LCL = Value - tmp * se$U,
                             UCL = Value + tmp * se$U,
                             Dataset = dataset_name,
                             Diversity = index) %>% select(c("Dataset", "Order.q", "Size", "Value", "s.e.", "LCL", "UCL", "Diversity"))
    
    V = V %>% mutate(        Size = n,
                             s.e. = se$V,
                             LCL = Value - tmp * se$V,
                             UCL = Value + tmp * se$V,
                             Dataset = dataset_name,
                             Diversity = index) %>% select(c("Dataset", "Order.q", "Size", "Value", "s.e.", "LCL", "UCL", "Diversity"))
    
    S = S %>% mutate(        Size = n,
                             s.e. = se$S,
                             LCL = Value - tmp * se$S,
                             UCL = Value + tmp * se$S,
                             Dataset = dataset_name,
                             Diversity = index) %>% select(c("Dataset", "Order.q", "Size", "Value", "s.e.", "LCL", "UCL", "Diversity"))
    
    
    if (diversity == "PD") {
      
      gamma = gamma %>% mutate(Reftime = reft)
      
      alpha = alpha %>% mutate(Reftime = reft)
      
      beta  = beta  %>% mutate(Reftime = reft)
      
      C     =  C    %>% mutate(Reftime = reft)
      
      U     =  U    %>% mutate(Reftime = reft)
      
      V     =  V    %>% mutate(Reftime = reft)
      
      S     =  S    %>% mutate(Reftime = reft)
    }
    
    if (diversity == "FD" & FDtype == "tau_value") {
      
      gamma = gamma %>% mutate(Tau = FDtau)
      
      alpha = alpha %>% mutate(Tau = FDtau)
      
      beta  = beta  %>% mutate(Tau = FDtau)
      
      C     =  C    %>% mutate(Tau = FDtau)
      
      U     =  U    %>% mutate(Tau = FDtau)
      
      V     =  V    %>% mutate(Tau = FDtau)
      
      S     =  S    %>% mutate(Tau = FDtau)
    }
    
    rbind(gamma, alpha, beta, C, U, V, S) %>% cbind(., Type = rep(c('gamma', 'alpha', 'beta', '1-C', '1-U', '1-V', '1-S'), each = length(q)))
    
  }
  
  output = lapply(1:length(data_list), function(i) for_each_dataset(data = data_list[[i]], dataset_name = dataset_names[i], N = Ns[i]))
  
  names(output) = dataset_names
  
  return(output)
  
}



#' ggplot for beta diversity
#' 
#' \code{ggObsbeta3D}: ggplot for observed beta diversity with order q
#' 
#' @param output the output from the function Obsbeta3D
#' @param type \code{type = 'B'} for plotting gamma, alpha, and beta diversity; \code{type = 'D'} for plotting four dissimilarity indices.
#' 
#' @return a figure for gamma, alpha, and Beta diversity, or a figure for four dissimilarity indices.
#' 
#' @examples
#' ## Taxonomic diversity for abundance data
#' data(beetle_abu)
#' output1 = Obsbeta3D(data = beetle_abu, diversity = 'TD', datatype = 'abundance', 
#'                     nboot = 10)
#' 
#' ggObsbeta3D(output1, type = 'B')
#' ggObsbeta3D(output1, type = 'D')
#' 
#' 
#' ## Phylogenetic Hill numbers for abundance data
#' data(beetle_abu)
#' data(beetle_tree)
#' output2 = Obsbeta3D(data = beetle_abu, diversity = 'PD', datatype = 'abundance', 
#'                     nboot = 10, conf = 0.95, 
#'                     PDtree = beetle_tree, PDreftime = NULL, PDtype = 'meanPD')
#' 
#' ggObsbeta3D(output2, type = 'B')
#' ggObsbeta3D(output2, type = 'D')
#' 
#' 
#' ## Functional diversity for abundance data under single threshold
#' data(beetle_abu)
#' data(beetle_distM)
#' output3 = Obsbeta3D(data = beetle_abu, diversity = 'FD', datatype = 'abundance', 
#'                     nboot = 10, conf = 0.95, 
#'                     FDdistM = beetle_distM, FDtype = 'tau_value', FDtau = NULL)
#' 
#' ggObsbeta3D(output3, type = 'B')
#' ggObsbeta3D(output3, type = 'D')
#' 
#' 
#' ## Functional diversity for abundance data with thresholds integrating from 0 to 1
#' data(beetle_abu)
#' data(beetle_distM)
#' output4 = Obsbeta3D(data = beetle_abu, diversity = 'FD', datatype = 'abundance', 
#'                     nboot = 10, conf = 0.95, 
#'                     FDdistM = beetle_distM, FDtype = 'AUC', FDcut_number = 30)
#' 
#' ggObsbeta3D(output4, type = 'B')
#' ggObsbeta3D(output4, type = 'D')
#' 
#' 
#' ## Taxonomic diversity for incidence data
#' data(beetle_inc)
#' output5 = Obsbeta3D(data = beetle_inc, diversity = 'TD', datatype = 'incidence_raw', 
#'                     nboot = 10, conf = 0.95)
#' 
#' ggObsbeta3D(output5, type = 'B')
#' ggObsbeta3D(output5, type = 'D')
#' 
#' 
#' ## Phylogenetic Hill numbers for incidence data
#' data(beetle_inc)
#' data(beetle_tree)
#' output6 = Obsbeta3D(data = beetle_inc, diversity = 'PD', datatype = 'incidence_raw', 
#'                     nboot = 10, conf = 0.95, 
#'                     PDtree = beetle_tree, PDreftime = NULL)
#' 
#' ggObsbeta3D(output6, type = 'B')
#' ggObsbeta3D(output6, type = 'D')
#' 
#' 
#' ## Functional diversity for incidence data under single threshold
#' data(beetle_inc)
#' data(beetle_distM)
#' output7 = Obsbeta3D(data = beetle_inc, diversity = 'FD', datatype = 'incidence_raw', 
#'                     nboot = 10, conf = 0.95, 
#'                     FDdistM = beetle_distM, FDtype = 'tau_value', FDtau = NULL)
#'
#' ggObsbeta3D(output7, type = 'B')
#' ggObsbeta3D(output7, type = 'D')
#' 
#' 
#' ## Functional diversity for incidence data with thresholds integrating from 0 to 1
#' data(beetle_inc)
#' data(beetle_distM)
#' output8 = Obsbeta3D(data = beetle_inc, diversity = 'FD', datatype = 'incidence_raw', 
#'                     nboot = 10, conf = 0.95, 
#'                     FDdistM = beetle_distM, FDtype = 'AUC', FDcut_number = 30)
#' 
#' ggObsbeta3D(output8, type = 'B')
#' ggObsbeta3D(output8, type = 'D')
#' 
#' @export
ggObsbeta3D = function(output, type = 'B'){
  
  if (type == 'B'){
    
    gamma = lapply(output, function(y) y %>% filter(Type == 'gamma')) %>% do.call(rbind,.) %>% mutate(Type = "Gamma") %>% as_tibble()
    alpha = lapply(output, function(y) y %>% filter(Type == 'alpha')) %>% do.call(rbind,.) %>% mutate(Type = "Alpha") %>% as_tibble()
    beta =  lapply(output, function(y) y %>% filter(Type == 'beta'))  %>% do.call(rbind,.) %>% mutate(Type = "Beta")  %>% as_tibble()
    
    df = rbind(gamma, alpha, beta)
    df$Type <- factor(df$Type, levels = c("Gamma","Alpha","Beta"))
    
  }
  
  if (type == 'D'){
    
    C = lapply(output, function(y) y %>% filter(Type == '1-C')) %>% do.call(rbind,.) %>% mutate(Type = "1-CqN") %>% as_tibble()
    U = lapply(output, function(y) y %>% filter(Type == '1-U')) %>% do.call(rbind,.) %>% mutate(Type = "1-UqN") %>% as_tibble()
    V = lapply(output, function(y) y %>% filter(Type == '1-V')) %>% do.call(rbind,.) %>% mutate(Type = "1-VqN") %>% as_tibble()
    S = lapply(output, function(y) y %>% filter(Type == '1-S')) %>% do.call(rbind,.) %>% mutate(Type = "1-SqN") %>% as_tibble()
    
    df = rbind(C, U, V, S)
    df$Type <- factor(df$Type, levels = c("1-CqN", "1-UqN", "1-VqN", "1-SqN"))
    
  }
  
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  
  if (unique(output[[1]]$Diversity) == 'TD'    ) { ylab = "Taxonomic diversity" }
  if (unique(output[[1]]$Diversity) == 'PD'    ) { ylab = "Phylogenetic diversity" }
  if (unique(output[[1]]$Diversity) == 'meanPD') { ylab = "Mean phylogenetic diversity" }
  if (unique(output[[1]]$Diversity) == 'FD_tau') { ylab = "Functional diversity (given tau)" }
  if (unique(output[[1]]$Diversity) == 'FD_AUC') { ylab = "Functional diversity (AUC)" }
  
  ggplot(data = df) +
    geom_line(aes(x = Order.q, y = Value, col = Dataset), size=1.2) + 
    geom_ribbon(aes(x = Order.q, ymin = LCL, ymax = UCL, fill = Dataset, col = NULL), alpha = 0.4) + 
    scale_colour_manual(values = cbPalette) + 
    scale_fill_manual(values = cbPalette) + 
    facet_grid(Type ~ ., scales = 'free_y') +
    theme_bw() + 
    theme(legend.position = "bottom", 
          legend.title = element_blank(),
          strip.text = element_text(size = 15, face = 'bold'),
          axis.title = element_text(hjust = 0.5, size = 15, face = 'bold'),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.box = "vertical",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-10, -10, -5, -10),
          legend.text = element_text(size = 13),
          plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
    labs(x = 'Order q', y = ylab) +
    guides(linetype = guide_legend(keywidth = 2.5))
}


bootstrap_population_multiple_assemblage = function(data, data_gamma, datatype){
  
  if (datatype == 'abundance'){
    
    S_obs = sum(data_gamma > 0)
    n = sum(data_gamma)
    f1 = sum(data_gamma == 1)
    f2 = sum(data_gamma == 2)
    f0_hat = ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 1)/n * f1^2/2/f2) %>% ceiling()
    
    output = apply(data, 2, function(x){
      
      p_i_hat = iNEXT.3D:::EstiBootComm.Ind(Spec = x)
      
      if(length(p_i_hat) != length(x)){
        
        p_i_hat_unobs = p_i_hat[(length(x)+1):length(p_i_hat)]
        p_i_hat_obs = p_i_hat[1:length(x)]
        p_i_hat = c(p_i_hat_obs, rep(0, f0_hat))
        candidate = which(p_i_hat==0)
        
        chosen = sample(x = candidate, size = min(length(p_i_hat_unobs), length(candidate)), replace = F)
        p_i_hat[chosen] = (1-sum(p_i_hat))/length(chosen)
        
        p_i_hat
        
      } else {
        
        p_i_hat = c(p_i_hat, rep(0, f0_hat))
        p_i_hat
        
      }
    })
    
  }
  
  if (datatype == 'incidence'){
    
    S_obs = sum(data_gamma > 0)
    t = data_gamma[1]
    Q1 = sum(data_gamma == 1)
    Q2 = sum(data_gamma == 2)
    Q0_hat = if ( Q2 == 0 ){( (t-1)/t ) * ( Q1*(Q1-1)/2 )} else {( (t-1)/t ) * ( (Q1^2) / (2*Q2) )} %>% ceiling
    
    output = apply(data, 2, function(x){
      
      pi_i_hat = iNEXT.3D:::EstiBootComm.Sam(Spec = x)
      
      if(length(pi_i_hat) != (length(x) - 1)){
        
        pi_i_hat_unobs = pi_i_hat[length(x):length(pi_i_hat)]
        pi_i_hat_obs = pi_i_hat[1:(length(x)-1)]
        pi_i_hat = c(pi_i_hat_obs, rep(0, Q0_hat))
        candidate = which(pi_i_hat == 0)
        
        chosen = sample(x = candidate, size = min(length(pi_i_hat_unobs), length(candidate)), replace = F)
        pi_i_hat[chosen] = unique(pi_i_hat_unobs)
        
        pi_i_hat
        
      } else {
        
        pi_i_hat = c(pi_i_hat, rep(0, Q0_hat))
        pi_i_hat
        
      }
    })
    
  }
  
  return(output)
  
}

Bootstrap_distance_matrix = function(data, distance_matrix, f0.hat, datatype){
  
  if (datatype == "incidence_freq") {
    n = data[1]
    X = data[-1]
    u = sum(data)
  } else if (datatype == "abundance") {
    n = sum(data)
    X = data
  }
  
  # n = sum(data)
  distance = as.matrix(distance_matrix)
  dij = distance
  # X = data
  
  F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
  F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])
  
  if (datatype == "abundance") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
  } else if (datatype == "incidence_freq") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-1)^2 * (F11^2)/(4* n* n* F22)), ((n-1)* (n-1)* (F11*(F11-0.01))/(4 *n * n)) )
  }
  
  if (f0.hat == 0) {
    d = dij
  } else if (f0.hat == 1) {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X)*f0.hat), length(X), f0.hat)
    
    d00 = matrix(0, f0.hat, f0.hat)
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  } else {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X)*f0.hat), length(X), f0.hat)
    
    fo.num = (f0.hat * (f0.hat-1) )/2
    d00 = matrix(0, f0.hat, f0.hat)
    d00[upper.tri(d00)] = (F00hat/2)/fo.num
    d00 <- pmax(d00, t(d00))###signmatrix
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  }
  
  return(d)
  
}

