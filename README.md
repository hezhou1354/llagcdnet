# R Package {llagcdnet}
## He Zhou (zhou1354@umn.edu)
LASSO, Elastic Net (Adaptive) and Folded Concave (SCAD) Penalized Least Squares, Logistic Regression, HHSVM, Squared Hinge SVM, Expectile Regression, Probit Regression and Composite Probit Regressino using a Fast (LLA-)GCD Algorithm.


To check do
```
cd package
R CMD build llagcdnet
R CMD check llagcdnet_*.tar.gz
```
The vignette may take some time, probably want
```
R CMD check llagcdnet_*.tar.gz --no-vignettes
```
(instead of without `--no-vignettes`) except for one last check before commit.

