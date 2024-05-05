#' Convert neuropsychological test scores to demographically adjusted norms.
#'
#' @param data a data frame containing the variables needed for the norming
#' process. The current version of the function does not accomodate missing
#' data. For best results, exclude cases with missing test scores or missing
#' demographics before applying this function.
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
#' included into normative formulas. Demographic variables should be numeric or
#' binary (0/1).
#' @param mfp.alpha a numeric value between 0 and 1 that sets significance level
#' for inclusion of demographic predictors into normative formula. Passed to the
#' \code{mfp2()} function (its \code{select} argument). Default value is 1 for
#' inclusion of all predictors regardless of their significance.
#' @param rnd.a a logical indicating whether the adjusted scores (T-scores)
#' should be rounded. Default is TRUE.
#' @param mean.a numeric value for the mean of adjusted score (T-score)
#' distribution. Default is 50.
#' @param sd.a numeric value for the standard deviation of adjusted score
#' (T-score) distribution. Default is 10.
#'
#' @details
#' The \code{score2adjust()} function can be used by neuropsychologists, who
#' wish to construct normative formulas for cognitive tests that adjust for
#' expected effects of demographic characteristics (e.g., age), using methods
#' described in Heaton et al. (2003 & 2009). The adjusted scores are sometimes
#' referred to as T-scores in the literature. The norming procedure makes use of
#' the \code{mfp2()} function from the \code{mfp2} package to explore nonlinear
#' associations between cognition and demographic variables. Detailed
#' description of the procedure are found in Umlauf et al. (2024). (Previous
#' versions of the function depended on \code{mfp} package.)
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
#' Umlauf A et al. (2024) Automated procedure for demographic adjustments on
#' cognitive test scores. <doi:10.1080/23279095.2023.2288231>
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
#' @importFrom stats as.formula
#' @importFrom stats predict
#' @importFrom stats qt
#' @importFrom stats residuals
#' @importFrom stats sd
#' @import mfp2
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
#'              demographics = c("age", "male"))
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
                        paste("fp2(", demographics[i],
                              ", df = 4, center = FALSE, select = mfp.alpha)",
                              sep = ""))
  }
  # define model formula
  RS.fmla <- paste(unlist(covs), collapse = " + ")
  fmla <- as.formula(paste(test.score, " ~ ", RS.fmla))

  model <- mfp2(fmla, data = data1[data1[, group.id] == control.id, ],
                family = "gaussian", select = mfp.alpha,
                maxits = 5, verbose = FALSE)
  predicted <- as.numeric(predict(model,
                                  newdata = data1[,c(test.score,demographics)]))
  s <- sd(residuals(model), na.rm = TRUE)
  z <- (data1[, test.score] - predicted)/s
  t <- z*sd.a + mean.a

  {
    if (rnd.a == FALSE) data1[, ts.name] <- t  # adjusted scores (T-scores)
    else data1[, ts.name] <- round(t)          # rounded adj. scores (integers)
  }

  # calculate 95% confidence intervals for the coefficients
  conflev <- qt(0.975, df = model$df.residual)
  coef.me <- conflev*sqrt(diag(summary(model)$cov.scaled))
  coef.lower <- (model$coefficients - coef.me)
  coef.upper <- (model$coefficients + coef.me)

  # Save transformations for the output
  trafo <- cbind(model$transformations, model$fp_terms)
  { if(!"power2" %in% names(trafo)) trafo$power2 <- NA }
  trafo$formula <- NA
  for(i in 1:nrow(trafo)) # this loop assumes that 'center = FALSE'
  {
    if(trafo[i, "selected"] == TRUE){
      x.i <- rownames(trafo)[i]
      if(trafo$shift[i] == 0) x.shift <- x.i
      else{
        shft.i <- ifelse(trafo$shift[i] < 0, as.character(trafo$shift[i]),
                         as.character(paste("+",trafo$shift[i], sep="")))
        x.shift <- paste("(",x.i,shft.i,")",sep = "")
      }

      if(trafo$scale[i] == 1) x.scale <- x.shift
      else x.scale <- paste(x.shift,"/",trafo$scale[i], sep="")

      if(is.na(trafo$power2[i]))
        fmla.i <- ifelse(trafo$power1[i] !=0,
                         paste("I((",x.scale,")^",trafo$power1[i],")", sep=""),
                         paste("I(log(",x.scale,"))", sep=""))
      else if(trafo$power1[i] != trafo$power2[i]){
        fmla.i.p1 <- ifelse(trafo$power1[i] !=0,
                            paste("I((",x.scale,")^",trafo$power1[i],")",
                                  sep=""),
                            paste("I(log(",x.scale,"))", sep=""))
        fmla.i.p2 <- ifelse(trafo$power2[i] !=0,
                            paste("I((",x.scale,")^",trafo$power2[i],")",
                                  sep=""),
                            paste("I(log(",x.scale,"))", sep=""))
        fmla.i <- paste(fmla.i.p1, fmla.i.p2, sep="+")
      }
      else if(trafo$power1[i] == 0 & trafo$power2[i] == 0){
        fmla.i.p1 <- paste("I(log(",x.scale,"))", sep="")
        fmla.i.p2 <- paste("I(log(",x.scale,")^2)", sep="")
        fmla.i <- paste(fmla.i.p1, fmla.i.p2, sep="+")
      }
      else if(trafo$power1[i] == trafo$power2[i]){
        fmla.i.p1 <- paste("I((",x.scale,")^",trafo$power1[i],")", sep="")
        fmla.i.p2 <- paste("I((",x.scale,")^",trafo$power2[i],
                           "log(",x.scale,"))", sep="")
        fmla.i <- paste(fmla.i.p1, fmla.i.p2, sep="+")
      }
      trafo$formula[i] <- fmla.i
    }
  }

  ### save norming info
  tab.out <- list(data1[, test.score], data1[, ts.name],
                  MFP.formulas = list(full.model.summary =
                                        cbind(summary(model)$coef,
                                              low95 = coef.lower,
                                              upp95 = coef.upper),
                                      transformations = trafo[,"formula",
                                                              drop=FALSE],
                                      coefficients = model$coefficients,
                                      residual.SD = s,
                                      powers = trafo[,c("power1", "power2")]))
  names(tab.out)[1:2] <- c(test.score, ts.name)

  #### Print output
  return(tab.out)
}  # end of score2adjust() function
