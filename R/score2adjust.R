#' Convert neuropsychological test scores to demographically adjusted norms.
#'
#' @param data a data frame containing the variables needed for the norming
#' process.
#' @param test.score a character string specifying the name of the test to be
#' normed, usually the output of the \code{raw2scaled()} function.
#' @param group.id a character string specifying the name of the variable
#' containing group identification (i.e. control vs exposed/test/risk). Ignored,
#' if all.controls = TRUE.
#' @param control.id a character string specifying the label of the control
#' group within group.id variable. Ignored, if all.controls = TRUE.
#' @param all.controls a logical indicating whether all observations should be
#' treated as controls. Overwrites group.id and control.id.
#' @param demographics a single or multiple character strings (concatenated by
#' \code{c()} function) specifying the names of demographic predictors to be
#' included into normative formulas.
#' @param mfp.alpha a numeric value between 0 and 1 that sets significance level
#' for inclusion of demographic predictors into normative formula. Passed to the
#' \code{mfp()} function (its \code{select} argument). Default value is 1 for
#' inclusion of all predictors regardless of their significance.
#' @param rnd.a a logical indicating whether the adjusted scores (T-scores)
#' should be rounded. Default is TRUE.
#' @param mean.a numeric value for the mean of adjusted score (T-score)
#' distribution.
#' @param sd.a numeric value for the standard deviation of adjusted score
#' (T-score) distribution.
#'
#' @details
#' The \code{score2adjust()} function can be used by neuropsychologists, who
#' wish to construct normative formulas for cognitive tests that adjust for
#' expected effects of demographic characteristics (e.g., age), using methods
#' described in Heaton et al. (2003 & 2009). The adjusted scores are sometimes
#' referred to as T-scores in the literature. The norming procedure makes use of
#' the \code{mfp()} function from the \code{mfp} package to explore nonlinear
#' associations between cognition and demographic variables. Detailed
#' description of the procedure will be found in Umlauf et al. (in revision).
#'
#' @return
#' A list consisting of 3 objects. The first two are vectors containing the
#' non-adjusted test scores and the calculated demographically adjusted scores.
#' The last item in the output list is also a list called \code{MFP.formulas}.
#' It contains the information for calculation of adjusted scores, including
#' variable transformations (if any), multiple fractional polynomial (MFP) model
#' coefficients, the standard deviation of residuals resulting from the MFP
#' modeling, and a matrix with number of rows equal to the number of predictors
#' and 2 columns containing powers (in numeric form) selected for variable
#' transformations.
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
#' @export score2adjust
#'
#' @examples
#' data(PsychTestData)
#' PsychTestData$scaledscore <- raw2scaled(data=PsychTestData, test="rawscore",
#'                                         test.min=0, test.max=36,
#'                                         test.better="High", group.id="group",
#'                                         control.id="control")[[2]]
#' score2adjust(data = PsychTestData, test.score = "scaledscore",
#'              group.id = "group", control.id = "control",
#'              demographics = c("age", "sex"))
score2adjust <- function (data = NULL, test.score = NULL,
                          group.id = NULL, control.id = NULL,
                          all.controls = FALSE, demographics = NULL,
                          mfp.alpha = 1, rnd.a = TRUE, mean.a = 50, sd.a = 10)
{
  #### Check availability of necessary data:
  ## Stop the program if necessary information is not provided
  if (is.null(data))
    stop("Missing 'data' argument. Please provide name of the data frame.")
  if (is.null(test.score))
    stop("Missing 'test.score' argument.
         Please provide name of the test score for norming.")

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

  #### TEST SCORES --> DEMOGRAPHICALLY ADJUSTED SCORE (e.g., T-score)
  #### MFP procedure, calculating T-scores
  ## 'select' is set to 1 to keep all given variables in the model
  ## 'maxits' is the maximum number of iterations for the backfitting stage
  ## 'rescale' uses re-scaling to show the parameters for covariates on their
  ##           original scale (want FALSE)
  ## 'verbose' will hide the selection process if FALSE, and show them if TRUE

  ###### Calculate norms
  ### Make a data frame for internal purposes
  data1 <- data[, c(group.id, demographics, test.score)]
  ### Create variable names for calculated adjusted scores
  ts.name <- "score_adj"       # name the adjusted score

  ### Run norming process
  # calculate T-scores
  n.covs <- length(demographics)  # count demographic variables to control for
  covs <- list()
  for(i in 1:n.covs)
  {
    covs[[i]] <- ifelse(is.factor(data[, demographics[i]]),
                        demographics[i],
                        paste("fp(", demographics[i], ", df = 4)", sep = ""))
  }
  # define model formula
  RS.fmla <- paste(unlist(covs), collapse = " + ")
  fmla <- as.formula(paste(test.score, " ~ ", RS.fmla))

  model <- mfp(fmla, data = data1[data1[, group.id] == control.id, ],
               family = gaussian, select = mfp.alpha,
               maxits = 5, rescale = FALSE, verbose = FALSE)
  predicted <- predict(model, newdata = data1)
  s <- sd(residuals(model), na.rm = TRUE)
  z <- (data1[, test.score] - predicted)/s
  t <- z*sd.a + mean.a

  {
    if (rnd.a == FALSE) data1[, ts.name] <- t  # adjusted scores (T-scores)
    else data1[, ts.name] <- round(t)          # rounded adj. scores (integers)
  }

  # calculate 95% confidence intervals for the coefficients
  conflev = qt(0.975, df = model$df.residual)
  coef.lower <- (model$coefficients - conflev*sqrt(diag(model$var)))
  coef.upper <- (model$coefficients + conflev*sqrt(diag(model$var)))

  ### save norming info
  tab.out <- list(data1[, test.score], data1[, ts.name],
                  MFP.formulas = list(full.model.summary =
                                        cbind(summary(model)$coef,
                                              low95 = coef.lower,
                                              upp95 = coef.upper),
                                      transformations = model$trafo,
                                      coefficients = model$coefficients,
                                      residual.SD = s,
                                      powers = model$powers))
  names(tab.out)[1:2] <- c(test.score, ts.name)

  #### Print output
  return(tab.out)
}  # end of score2adjust() function
