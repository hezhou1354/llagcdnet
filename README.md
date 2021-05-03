# STAT 8054 Package {llagcdnet}
## He Zhou (zhou1354@umn.edu)



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

