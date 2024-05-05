## Updates to test2norm 0.3.0

* The new version of the package now depends on mfp2 package. Previous versions
  used mfp package.
* Functions test2norm() and score2adjust() were edited to use mfp2() and fp2()
  functions from the mfp2 package. 
* Functions test2norm() and score2adjust() may results in error when there are 
  missing test scores or demographic variables. Users are strongly advised to
  exclude cases with missing data prior to applying these functions.
* This version of the package accepts numeric and binary (coded as 0/1) 
  demographic predictors.

## Updates to test2norm 0.2.0

* Added a `NEWS.md` file to track changes to the package.
* Two functions were added: raw2scaled() and score2adjust(). These functions
separately perform 2 steps that are done within the test2norm() function.
Function raw2scaled() converts raw test scores to scaled scores and function
score2adjust() calculates demographically adjusted scores. 


