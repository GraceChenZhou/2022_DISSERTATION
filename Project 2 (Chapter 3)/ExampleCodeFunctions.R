################################################################################
# Title: Multilevel Bayesian joint modeling for a longitudinal marker and time-to-recurrent
# event to characterize heterogeneity in multi-center studies
# Copyright: Grace C. Zhou, Seongho Song, Pedro M. Afonso, Eleni-Rosalina Andrinopoulou, Rhonda D. Szczesniak
# Credit: Samuel Brilleman
# Date: 2022.06.30
# Notes: Define functions for standata.list
################################################################################

# Return a list (or vector if unlist = TRUE) which
# contains the embedded elements in list x named y 
fetch <- function(x, y, z = NULL, zz = NULL, null_to_zero = FALSE, 
                  pad_length = NULL, unlist = FALSE) {
  ret <- lapply(x, `[[`, y)
  if (!is.null(z))
    ret <- lapply(ret, `[[`, z)
  if (!is.null(zz))
    ret <- lapply(ret, `[[`, zz)
  if (null_to_zero) 
    ret <- lapply(ret, function(i) ifelse(is.null(i), 0L, i))
  if (!is.null(pad_length)) {
    padding <- rep(list(0L), pad_length - length(ret))
    ret <- c(ret, padding)
  }
  if (unlist) unlist(ret) else ret
}

# Drop intercept from a vector/matrix/array of named coefficients
drop_intercept <- function(x) { 
  
  # Check whether a vector/matrix/array contains an "(Intercept)"
  check_for_intercept <- function(x, logical = FALSE) {
    nms <- if (is.matrix(x)) colnames(x) else names(x)
    sel <- which("(Intercept)" %in% nms)
    if (logical) as.logical(length(sel)) else sel
  }
  
  sel <- check_for_intercept(x)
  if (length(sel) && is.matrix(x)) {
    x[, -sel, drop = FALSE]
  } else if (length(sel)) {
    x[-sel]
  } else {
    x
  }
}

# Create a named list using specified names or, if names are omitted, using the
# names of the objects in the list
nlist <- function(...) {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names) FALSE else nzchar(names(out))
  if (all(has_name)) 
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  } 
  
  return(out)
}

# Split the random effects part of a model formula into
#   - the formula part (ie. the formula on the LHS of "|"), and 
#   - the name of the grouping factor (ie. the variable on the RHS of "|")

split_at_bars <- function(x) {
  terms <- strsplit(deparse(x, 500), "\\s\\|\\s")[[1L]]
  if (!length(terms) == 2L)
    stop2("Could not parse the random effects formula.")
  re_form <- formula(paste("~", terms[[1L]]))
  group_var <- terms[[2L]]
  nlist(re_form, group_var)
}

# Return design matrices for the longitudinal submodel
make_assoc_parts_for_stan <- function(newdata,formulaLong,include_Zt = FALSE) {
  
  # construct model frame using predvars
  ##re design matrices
  model_frame <- stats::model.frame(lme4::subbars(formulaLong),newdata)
  x_form <- lme4::nobars(formulaLong)
  x <- model.matrix(x_form,model_frame)
  predictors.x <- drop_intercept(x)
  x.bar <- colMeans(predictors.x)
  x.scl <- sweep(predictors.x,2,x.bar, FUN='-')
  
  ## fe design matrices
  bars <- lme4::findbars(formulaLong)
  z_parts <- lapply(bars, split_at_bars)
  z_forms <- fetch(z_parts, "re_form")
  z <- lapply(z_forms, model.matrix, model_frame)
  names(z) <- fetch(z_parts,'group_var')
  ret <- nlist(x, predictors.x, x.scl, x.bar, z) # return list
  
  # optionally add the sparse Zt matrix
  if (include_Zt) 
    ret$Zt <- lme4::mkReTrms(bars, model_frame)$Zt
  
  ret
}

# Function to return standardised GK quadrature points and weights
get_quadpoints <- function(nodes = 15) {
  if (!is.numeric(nodes) || (length(nodes) > 1L)) {
    stop("'qnodes' should be a numeric vector of length 1.")
  } else if (nodes == 15) {
    list(
      points = c(
        -0.991455371120812639207,
        -0.949107912342758524526,
        -0.86486442335976907279,
        -0.7415311855993944398639,
        -0.5860872354676911302941,
        -0.4058451513773971669066,
        -0.2077849550078984676007,
        0,
        0.2077849550078984676007,
        0.405845151377397166907,
        0.5860872354676911302941,
        0.741531185599394439864,
        0.86486442335976907279,
        0.9491079123427585245262,
        0.991455371120812639207),
      weights = c(
        0.0229353220105292249637,
        0.063092092629978553291,
        0.10479001032225018384,
        0.140653259715525918745,
        0.1690047266392679028266,
        0.1903505780647854099133,
        0.204432940075298892414,
        0.209482141084727828013,
        0.204432940075298892414,
        0.1903505780647854099133,
        0.169004726639267902827,
        0.140653259715525918745,
        0.1047900103222501838399,
        0.063092092629978553291,
        0.0229353220105292249637))      
  } else if (nodes == 11) {
    list(
      points = c(
        -0.984085360094842464496,
        -0.906179845938663992798,
        -0.754166726570849220441,
        -0.5384693101056830910363,
        -0.2796304131617831934135,
        0,
        0.2796304131617831934135,
        0.5384693101056830910363,
        0.754166726570849220441,
        0.906179845938663992798,
        0.984085360094842464496),
      weights = c(
        0.042582036751081832865,
        0.1152333166224733940246,
        0.186800796556492657468,
        0.2410403392286475866999,
        0.272849801912558922341,
        0.2829874178574912132043,
        0.272849801912558922341,
        0.241040339228647586701,
        0.186800796556492657467,
        0.115233316622473394025,
        0.042582036751081832865))     
  } else if (nodes == 7) {
    list(
      points = c(
        -0.9604912687080202834235,
        -0.7745966692414833770359,
        -0.4342437493468025580021,
        0,
        0.4342437493468025580021,
        0.7745966692414833770359,
        0.9604912687080202834235),
      weights = c(
        0.1046562260264672651938,
        0.268488089868333440729,
        0.401397414775962222905,
        0.450916538658474142345,
        0.401397414775962222905,
        0.268488089868333440729,
        0.104656226026467265194))      
  } else stop("'qnodes' must be either 7, 11 or 15.")  
}

# Unlist the result from an lapply call
uapply <- function(X, FUN, ...) {
  unlist(lapply(X, FUN, ...))
}

# Convert a standardised quadrature node to an unstandardised value based on 
# the specified integral limits
unstandardise_qpts <- function(t, a, b) {
  ((b - a) / 2) * t + ((b + a) / 2)
}


# Convert a standardised quadrature weight to an unstandardised value based on 
# the specified integral limits
unstandardise_qwts <- function(t, a, b) {
  ((b - a) / 2) * t
}

# Construct a list with information on the event submodel
handle_e_mod <- function(data, qnodes, id_var, width) {
  
  qwts0 <- list()
  qpts0.gap <- list()
  qpts0.cal <- list()
  
  qq <- get_quadpoints(qnodes)
  
  for (i in unique(data$id)){
    
    subdata=data %>% filter(id==i)
    Nobs=nrow(subdata)
    
    for (j in 1:Nobs){
      
      startime=subdata$tstart[j]
      stoptime=subdata$tstop[j]
      
      qwts0[[paste(i,j,sep='_')]] <- uapply(qq$weights, unstandardise_qwts, startime, stoptime)
      qpts0.gap[[paste(i,j,sep='_')]] <- uapply(qq$points, unstandardise_qpts, 0, stoptime-startime) #gap time scale
      qpts0.cal[[paste(i,j,sep='_')]] <- uapply(qq$points, unstandardise_qpts, startime, stoptime) #calendar time scale
      
    }
  }
  
  remove(startime,stoptime,i,j)
  
  qids <- rep(data$id,each=qnodes)
  qwts <- unlist(qwts0)
  qpts.gap <- unlist(qpts0.gap)
  qpts.cal <- unlist(qpts0.cal)
  
  remove(qwts0,qpts0.gap,qpts0.cal)
  
  # Event times/ids (for failures only)
  epts.gap <- data %>% filter(status==1) %>% mutate(gap.t=tstop-tstart) %>% pull(gap.t)# event times (for failures only)
  epts.cal <- data %>% filter(status==1) %>% pull(tstop)# event times (for failures only)
  eids <- data %>% filter(status==1) %>% pull(id) # subject ids (for failures only)
  
  # Both event times/ids and quadrature times/ids
  cpts.gap <- c(epts.gap, qpts.gap)
  cpts.cal <- c(epts.cal, qpts.cal)
  cids <- c(eids, qids) # NB using c(.) demotes factors to integers
  
  # Mean log incidence rate - used for shifting log baseline hazard
  norm_const <- log(sum(data$status) / sum(data$tstop))
  
  nlist(Npat = length(unique(data$id)), Nevents = sum(data$status), 
        qnodes, qwts, qpts.gap, qpts.cal, 
        qids, epts.gap, epts.cal, eids, cpts.gap, 
        cpts.cal, cids, norm_const)
  
}

# Function to return the range or SD of the predictors, used for scaling the priors
get_scale_value <- function(x) {
  num.categories <- n_distinct(x)
  x.scale <- 1
  if (num.categories == 2) {
    x.scale <- diff(range(x))
  } else if (num.categories > 2) {
    x.scale <- sd(x)
  }
  return(x.scale)
}

