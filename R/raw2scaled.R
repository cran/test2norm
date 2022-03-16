#' Convert raw neuropsychological test scores to scaled scores.
#'
#' @param data a data frame containing the test score
#' @param test a character string specifying the name of the variable containing
#' the test score
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
#' @param rnd.s a logical indicating whether the scaled scores should be
#' rounded. Default is FALSE.
#'
#' @details
#' The \code{raw2scaled()} function can be used by neuropsychologists, who wish
#' to convert raw test scores to scaled scores (mean=10, SD=3), using methods
#' described in Heaton et al. (2003 & 2009). The raw test scores that have many
#' decimal digits should be rounded to fewer digits prior to the application of
#' the \code{raw2scaled()} function. This will significantly reduce software
#' running time. The recommended number of decimal digits is 4 or fewer. Values
#' below \code{test.min} or above \code{test.max} will result in NA. Detailed
#' description of the procedure will be found in Umlauf et al. (In revision).
#'
#' Note that the function does not guarantee that the calculated scaled scores
#' range from 0 to 20, because it uses normal distribution quantiles. Therefore,
#' scaled scores outside 0-20 range are possible.
#'
#' @return
#' A list consisting of 3 objects. The first two are vectors containing the
#' original raw test scores and the calculated scaled scores. The third object
#' in the list, called \code{SS.maps}, contains conversions from raw scores to
#' scaled scores in a form of a table with two columns, one representing scaled
#' scores (one per row) and one representing raw scores (a single value or range
#' of raw values corresponding to each scaled score). Note that this table
#' shows rounded scaled scores regardless of the value for \code{rnd.s}.
#'
#' @author
#' Anya Umlauf
#'
#' @references
#' Umlauf A et al. (2022) Automated procedure for demographic adjustments on
#' cognitive test scores. Manuscript submitted for publication.
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
#' @import stats
#' @export raw2scaled
#'
#' @examples
#' data(PsychTestData)
#' raw2scaled(data = PsychTestData, test = "rawscore",
#'            test.min = 0, test.max = 36, test.better = "High",
#'            group.id = "group", control.id = "control")
raw2scaled <- function (data = NULL, test = NULL, test.min = NULL,
                        test.max = NULL, test.better = c("High", "Low"),
                        group.id = NULL, control.id = NULL,
                        all.controls = FALSE, rnd.s = FALSE)
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
  { data$group <- "control"; group.id <- "group"; control.id <- "control" }
  else
  {
    if (is.null(group.id))
      stop("Missing 'group.id' argument. User must provide grouping variable.")
    if (!is.null(group.id) & is.null(control.id))
      stop("Missing 'control.id' value.")
  }

  #### function to count number of decimal places, needed for SS mappings
  decimals <- function(x)
  {
    k <- 0
    # search for largest number of decimal places seen in the data
    while(!isTRUE(all.equal(x, round(x, k)))) k <- k+1
    return(k)
  }

  #### RAW TEST SCORES --> SCALED SCORES
  #### Scaled score procedure
  ## 'raw' is a vector of raw test scores for all subjects
  ## 'control' is a vector of raw test scores for controls
  ## 'direction' = High if larger values of a test mean better results
  ## 'direction' = Low if larger values of a test mean worse results
  ## 'direction' is needed to achieve SS always scaled from worse to better
  ## 'minRaw' & 'maxRaw' are possible raw score limits (for consistency
  ## between cohorts and for a 'nicer' distribution at the ends
  ## of the scaled scores (to use extrapolation))
  scale.scores = function(raw, control, direction, minRaw, maxRaw)
  {
    ## step 1: calculate scaled scores (SS) for controls
    # use only those controls who have scores in the plausible range
    controls <- control[!is.na(control) & control <= maxRaw & control >= minRaw]
    n <- length(controls)
    ranks <- rank(controls, ties.method = "average")
    p.ties <- ranks/(n+1)

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

  ###### Calculate scaled scores, applying above functions
  ### Make a data.frame for internal purposes
  data1 <- data[, c(group.id, test)]
  ### Create a name for calculated score
  ss.name <- paste(test, "scaled", sep = "_")     # name the scaled score

  ### Calculate scaled scores
  ss <- scale.scores(raw = data[, test],
                     control = data[data[, group.id] == control.id, test],
                     direction = test.better,
                     minRaw = test.min, maxRaw = test.max)
  # save non-rounded or rounded (integer) scaled scores
  {
    if (rnd.s == FALSE) data1[, ss.name] <- ss[[1]]  # SSs
    else data1[, ss.name] <- round(ss[[1]])          # rounded SSs (integers)
  }

  ### save scaled score info
  tab.out <- list(data1[, test], data1[, ss.name], SS.maps = ss[[2]])
  names(tab.out)[1:2] <- c(test, ss.name)

  #### Print output
  return(tab.out)
}  # end of raw2scaled() function
