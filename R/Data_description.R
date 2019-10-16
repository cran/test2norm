#' Neuropsychological test data
#'
#' A simulated data containing raw test scores and demographic characteristics
#' for 250 persons, 200 in the control group and 50 in the test group.
#' The raw test scores are to be converted to demographically corrected normed
#' scores, adjusting for effects of age and sex. The control group is used to
#' generate the norming formulas, which are then applied to all scores.

#'
#' @format A data frame with 250 rows and 4 variables:
#' \describe{
#' \item{rawscore}{raw test score on a neuropsychological test, ranging 0-36,
#' with higher values indicating better test performance}
#' \item{age}{age of the participant, in years}
#' \item{sex}{sex of the participant, male (M) or female (F)}
#' \item{group}{norming group the participant belongs to (control or test)}
#' }
#' @examples
#' data(PsychTestData)
#' test2norm(data = PsychTestData, test = "rawscore",
#'           test.min = 0, test.max = 36, test.better = "High",
#'           group.id = "group", control.id = "control",
#'           demographics = c("age", "sex"))
"PsychTestData"


