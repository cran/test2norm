#' Convert raw neuropsychological test scores to demographically adjusted norms.
#'
#' @param data a data frame containing the variables needed for the norming
#' process
#' @param test a character string specifying the name of the test to be normed
#' @param test.min a real number indicating the smallest possible test score
#' @param test.max a real number indicating the largest possible test score
#' @param test.better a character string indicating direction of the scores.
#' Use "High" if high test scores imply better performance, use "Low" otherwise.
#' @param group.id a character string specifying the name of the variable
#' containing group identification (i.e. control vs exposed/test/risk). Ignored,
#' if all.controls = TRUE.
#' @param control.id a character string specifying the label of the control
#' group within group.id variable. Ignored, if all.controls = TRUE.
#' @param all.controls a logical indicating whether all observations should be
#' treated as controls. Overwrites group.id and control.id.
#' @param demographics a single or multiple character strings (concatenated by
#' c() function) specifying the names of demographic predictors to be included
#' into normative formulas.
#' @param mfp.alpha a numeric value between 0 and 1 that sets significance level
#' for inclusion of demographic predictors into normative formula. Passed to the
#' mfp() function. Default value is 1 for inclusion of all predictors,
#' regardless of their significance.
#' @param rnd.s a logical indicating whether the scaled scores should be
#' rounded. Default is TRUE.
#' @param rnd.a a logical indicating whether the adjusted scores (T-scores)
#' should be rounded. Default is TRUE.
#' @param mean.a numeric value for the mean of adjusted score (T-score)
#' distribution.
#' @param sd.a numeric value for the standard deviation of adjusted score
#' (T-score) distribution.
#'
#' @details
#' The \code{test2norm()} function can be used by neuropsychologists, who wish
#' to construct normative formulas for cognitive tests that adjust for expected
#' effects of demographic characteristics (e.g., age), using methods described
#' in Heaton et al. (2003 & 2009). The norming procedure makes use of the
#' \code{mfp()} function from the \code{mfp} package to explore nonlinear
#' associations between cognition and demographic variables. The raw test scores
#' that have many decimal digits should be rounded to fewer digits prior to the
#' application of the \code{test2norm()} function. This will significantly
#' reduce software running time. The recommended number of decimal digits is 4
#' or fewer. Detailed description of the procedure will be found in Umlauf et al
#' (2019).
#'
#' @return
#' A list consisting of 6 objects. The first four are vectors containing the
#' original raw test scores and the calculated scaled scores, demographically
#' adjusted scores, and deficit scores. The fifth object in the list, called
#' \code{SS.maps}, contains conversions from raw scores to scaled scores in a
#' form of a table with two columns, one representing scaled scores (one per
#' row) and one representing raw scores (range of raw values corresponding to
#' each scaled score). The last item in the output list is also a list called
#' \code{MFP.formulas} and contains the information for calculation of adjusted
#' scores, including variable transformations (if any), multiple fractional
#' polynomial (MFP) model coefficients, and the standard deviation of residuals
#' resulting from the MFP modeling.
#'
#' @author
#' Anya Umlauf
#'
#' @references
#' Umlauf A et al (2019) Automated procedure to produce normative correction
#' formulas modeling demographic effects on cognitive test scores and apply them
#' to obtain demographically corrected scores. Manuscript submitted for
#' publication.
#'
#' Heaton RK, Taylor MJ, & Manly J (2003) Demographic effects and use of
#' demographically corrected norms with the WAIS-III and WMS-III. In: Tulsky D
#' et al. (Eds.) \emph{Clinical Interpretation of the WAIS-III and WMS-III}.
#' San Diego, CA: Academic Press, 183-210.
#'
#' Heaton RK, Ryan L, & Grant I (2009) Demographic influences and use of
#' demographically corrected norms in neuropsychological assessment. In Grant I
#' & Adams KM (Eds.) \emph{Neuropsychological Assessment of Neuropsychiatric
#' and Neuromedical Disorders}. New York, NY: Oxford University Press, 127-155.
#'
#' Benner A (2005) mfp: Multivariable fractional polynomials.
#' \emph{R News} 5(2): 20–23.
#'
#' @import stats
#' @import mfp
#' @export test2norm
#'
#' @examples
#' data(PsychTestData)
#' test2norm(data = PsychTestData, test = "rawscore",
#'           test.min = 0, test.max = 36, test.better = "High",
#'           group.id = "group", control.id = "control",
#'           demographics = c("age", "sex"))
test2norm <- function (data = NULL, test = NULL, test.min = NULL,
                       test.max = NULL, test.better = c("High", "Low"),
                       group.id = NULL, control.id = NULL,
                       all.controls = FALSE, demographics = NULL, mfp.alpha = 1,
                       rnd.s = TRUE, rnd.a = TRUE, mean.a = 50, sd.a = 10)
{
  #### Check availability of necessary data:
  ## Stop the program if necessary information is not provided
  if (is.null(data))
    stop("Missing 'data' argument. Please provide name of the data frame.")
  if (is.null(test))
    stop("Missing 'test' argument.
         Please provide name of the test for norming.")
  if (is.null(test.min))
    stop("Missing 'test.min' argument.
         Please provide minimum possible value for the test.")
  if (is.null(test.max))
    stop("Missing 'test.max' argument.
         Please provide maximum possible value for the test.")
  if (!is.numeric(test.min))
    stop("'test.min' is missing or not numeric")
  if (!is.numeric(test.max))
    stop("'test.max' in test info file is missing or not numeric")
  if (is.null(test.better))
    stop("Missing 'test.better' argument. Please specify 'High' or 'Low'.")
  if (!test.better %in% c("Low", "High"))
    stop("'test.better' should be equal to 'Low' or 'High'")

  if (is.logical (all.controls) & all.controls == TRUE)
    { data$group <- "control"; group.id <- "group"; control.id <- "control"}
  else
  {
	  if (is.null(group.id))
	    stop("Missing 'group.id' argument. User must provide grouping variable.")
	  if (!is.null(group.id) & is.null(control.id))
	    stop("Missing 'control.id' value.")
  }
  if (is.null(demographics)) stop("Missing 'demographics' argument.")

  #### function used to count number of decimal places, needed for SS mappings
  decimals <- function(x)
  {
    k <- 0
    # search for largest number of decimal places seen in the data
    while(!isTRUE(all.equal(x, round(x, k)))) k <- k+1
    return(k)
  }

  # The norms calculation is executed in three steps:
  # STEP 1: raw test scores --> scaled scores
  # STEP 2: scaled scores --> demographically corrected T-scores
  # STEP 3: T-scores --> deficit scores

  #### STEP 1: raw test scores --> scaled scores
  #### Scale score procedure
    ## 'raw' is a vector of raw test scores for all subjects
		## 'control' is a vector of raw test scores for controls
		## 'direction' = High if larger values of a test mean better results
		## 'direction' = Low if larger values of a test mean worse results
		## 'direction' is needed to achieve SS always scaled from worse to better
		## 'minRaw' & 'maxRaw' are possible raw score limits (for consistency
		## between cohorts and for a 'nicer' distribution at the ends
    ## of the scalled scores (to use extrapolation)
  scale.scores = function(raw, control, direction, minRaw, maxRaw)
  {
    ## step 1: calculate scaled scores (SS) for controls
    # use only those controls who have scores in the plausible range
  	controls <- control[!is.na(control) & control <= maxRaw & control >= minRaw]
    n <- length(controls)
    ranks <- rank(controls, ties.method = "average")
	  p.ties <- ranks/(n+1) # normal quantiles

	  if (direction == "Low") quantiles <- qnorm(1 - p.ties, mean = 0, sd = 1)
	  else quantiles <- qnorm(p.ties, mean = 0, sd = 1)

    scale.contr <- scale(quantiles)*3+10
    meanQ <- mean(quantiles, na.rm = TRUE)
    sdQ <- sd(quantiles, na.rm = TRUE)

    if (min(controls, na.rm = TRUE) != minRaw)
    { controls <- c(controls, minRaw)
      if (direction == "Low")
        quantileMIN <- qnorm((1-0.5/(n+1)), mean = 0, sd = 1)
      else quantileMIN <- qnorm((0.5/(n+1)), mean = 0, sd = 1)
      scale.contr <- c(scale.contr, 10 + 3*(quantileMIN - meanQ)/sdQ)
    }
    if (max(controls, na.rm = TRUE) != maxRaw)
    { controls <- c(controls, maxRaw)
      if (direction == "Low")
        quantileMAX <- qnorm((0.5/(n+1)), mean = 0, sd = 1)
      else quantileMAX <- qnorm((1-0.5/(n+1)), mean = 0, sd = 1)
      scale.contr <- c(scale.contr, 10 + 3*(quantileMAX - meanQ)/sdQ)
    }

  ## step 2: use SS from controls to calculate SS for the all
    scaled <- c()
    cont <- data.frame(rs = controls, ss = scale.contr)
    cont <- na.omit(unique(cont[order(cont$rs), ]))

    raw[!is.na(raw) & raw > maxRaw] <- NA # exclude values outside of range
    raw[!is.na(raw) & raw < minRaw] <- NA # exclude values outside of range

    scaled <- scale.contr[match(raw, controls)]

    for(i in 1:length(raw)) # calculate SSs for values not observed in controls
    {
    	if (is.na(scaled[i]) && !is.na(raw[i]))
        {
        	k <- length(cont$rs[raw[i] > cont$rs])               # position
        	w <- (cont$rs[k+1]-raw[i])/(cont$rs[k+1]-cont$rs[k]) # weight
        	scaled[i] <- cont$ss[k]*w + cont$ss[k+1]*(1-w)       # weighted SS
        }
    }

  ## step 3: use SS from controls to create mapping files for each test
    ndp <- decimals(round(controls, 4))			   # ndp (number of decimal places)
    raw.map <- seq(minRaw, maxRaw, by = 1/(10^ndp))
    scale.map <- scale.contr[match(raw.map, controls)]

    for(i in 1:length(scale.map))
    {
    	if (is.na(scale.map[i]))
        {
        	k <- length(cont$rs[raw.map[i] > cont$rs])
        	w <- (cont$rs[k+1]-raw.map[i])/(cont$rs[k+1]-cont$rs[k])
        	scale.map[i] <- cont$ss[k]*w + cont$ss[k+1]*(1-w)
        }
    }

	scale.map.round <- round(scale.map)       	# a vector of rounded SSs
	Uscale.map.round <- unique(scale.map.round)	# a vector of unique rounded SSs
	mapping <- data.frame(raw = NA, ss = NA)  	# mapping table

	for(i in 1:length(Uscale.map.round))
	{
		interval <- raw.map[scale.map.round == Uscale.map.round[i]]
		mapping[i,] <- c(ifelse(length(interval) == 1,
		                        formatC(min(interval), digits = ndp, format = "f"),
		                        paste(formatC(min(interval),
		                                      digits = ndp, format = "f"), " - ",
		                              formatC(max(interval),
		                                      digits = ndp, format = "f"))),
		                 Uscale.map.round[i])
	}

	sc.scores <- list(scaled, mapping)
  return(sc.scores)

}   ## END of scaled score calculating function

#### STEP 2: at the moment is done within 'for' loop down below
#### MFP procedure, calculating T-scores
  ## 'select' is set to 1 to keep all given variables in the model
	## 'maxits' is the maximum number of iterations for the backfitting stage
	## 'rescale' uses re-scaling to show the parameters
  ## for covariates on their original scale (want FALSE)
	## 'verbose' will hide the selection process if FALSE, and show them if TRUE
## END of t-score calculating function

#### STEP 3:
#### Deficit scores procedure
deficit.scores <- function(Tscore)
{
	def <- rep(NA, length(Tscore))
	for (i in 1:length(Tscore))
	{
	  if (!is.na(Tscore[i]))
	  {
	    if (Tscore[i] >= (mean.a - sd.a)) {def[i] <- 0}
	    else if (Tscore[i] >= (mean.a - 1.5*sd.a)) {def[i] <- 1}
	    else if (Tscore[i] >= (mean.a - 2.0*sd.a)) {def[i] <- 2}
	    else if (Tscore[i] >= (mean.a - 2.5*sd.a)) {def[i] <- 3}
	    else if (Tscore[i] >= (mean.a - 3.0*sd.a)) {def[i] <- 4}
	    else if (Tscore[i] < (mean.a - 3.0*sd.a))  {def[i] <- 5}
	  }
	}
	return(def)

}## END of deficit score calculating function

  ###### Calculate norms, applying above functions
  ### Make a data.frame for internal purposes
  data1 <- data[, c(group.id, demographics, test)]
  ### Count observations
  Ntotal <- nrow(data1) # count sample size (# of unique subjects)
  Ncontrol <- nrow(data1[data1[,group.id] == control.id, ]) 	# count controls
  ### create a list to collect all necessary info
  ss.name <- paste(test, "scaled", sep = "_")     # name the scaled score
  ts.name <- paste(test, "adj", sep = "_")       # name the t-score
  tsr.name <- paste(test, "adj_rnd", sep = "_")  # name the rounded t-score
  ds.name <- paste(test, "def", sep = "_")        # name the deficit score

  ### run norming process
  # Step 1: calculate scaled scores
  info <- data.frame(test, test.min, test.max, test.better) ## test information
  ss <- scale.scores(raw = data[, test],
                     control = data[data[, group.id] == control.id, test],
                     direction = info$test.better,
                     minRaw = info$test.min, maxRaw = info$test.max)
  {
	  if (rnd.s == FALSE) data1[, ss.name] <- ss[[1]]  # SSs
	  else data1[, ss.name] <- round(ss[[1]])          # rounded SSs (integers)
  }

  # Step 2: calculate T-scores
  n.covs <- length(demographics)  # count demographic variables to control for
  covs <- list()
  for(i in 1:n.covs)
  {
	  covs[[i]] <- ifelse(is.factor(data[, demographics[i]]),
	                      demographics[i],
	                      paste("fp(", demographics[i], ", df = 4)", sep = ""))
  }
  RS.fmla <- paste(unlist(covs), collapse = " + ")
  fmla <- as.formula(paste(ss.name, " ~ ", RS.fmla))

  model <- mfp(fmla, data = data1[data1[, group.id] == control.id, ],
               family = gaussian, select = mfp.alpha,
               maxits = 5, rescale = FALSE, verbose = FALSE)
  predicted <- predict(model, newdata = data1)
  s <- sd(residuals(model), na.rm = TRUE)
  z <- (data1[, ss.name] - predicted)/s
  t <- z*sd.a + mean.a

  {
	  if (rnd.a == FALSE) data1[, ts.name] <- t  # T-scores
	  else data1[, ts.name] <- round(t)          # rounded T-scores (integers)
  }

  # Step 3: calculate Deficit Scores (DS)
  data1[, ds.name] <- deficit.scores(data1[, ts.name]) # DSs

  # Additional step: calculate 95% confidence intervals for the coefficients
  conflev = qt(0.975, df = model$df.residual)
  coef.lower <- (model$coefficients - conflev*sqrt(diag(model$var)))
  coef.upper <- (model$coefficients + conflev*sqrt(diag(model$var)))

  ### save norming info
  tab.out <- list(data1[, test], data1[, ss.name],
                  data1[, ts.name], data1[, ds.name],
                  SS.maps = ss[[2]],
                  MFP.formulas = list(full.model.summary =
                                        cbind(summary(model)$coef,
                                              low95 = coef.lower,
                                              upp95 = coef.upper),
                                      transformations = model$trafo,
                                      coefficients = model$coefficients,
                                      residual.SD = s))
  names(tab.out)[1:4] <- c(test, ss.name, ts.name, ds.name)

  #### Print output
  return(tab.out)
}  # end of test2norm() function
