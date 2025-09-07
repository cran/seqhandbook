# seqhandbook

## Miscellaneous Tools for Sequence Analysis

<!-- badges: start -->
[![R-CMD-check](https://framagit.org/nicolas-robette/seqhandbook/badges/master/pipeline.svg?key_text=R+CMD+check&key_width=90)](https://framagit.org/nicolas-robette/seqhandbook/-/commits/master)
[![](https://img.shields.io/gitlab/last-commit/nicolas-robette%2Fseqhandbook?gitlab_url=https%3A%2F%2Fframagit.org)](https://img.shields.io/gitlab/last-commit/nicolas-robette%2Fseqhandbook?gitlab_url=https%3A%2F%2Fframagit.org)
[![](https://www.r-pkg.org/badges/version/seqhandbook?color=blue)](https://cran.r-project.org/package=descriptio)
[![](https://www.r-pkg.org/badges/last-release/seqhandbook?color=blue)](https://cran.r-project.org/package=seqhandbook)
[![](https://img.shields.io/badge/DOI-10.32614/CRAN.package.seqhandbook-1f57b6?style=flat&link=https://doi.org/10.32614/CRAN.package.seqhandbook)](https://doi.org/10.32614/CRAN.package.seqhandbook)[![](http://cranlogs.r-pkg.org/badges/last-month/seqhandbook?color=orange)](https://cran.r-project.org/package=seqhandbook)
[![](http://cranlogs.r-pkg.org/badges/grand-total/seqhandbook?color=orange)](https://cran.r-project.org/package=seqhandbook)
<!-- badges: end -->


This R package complements the handbook on sequence analysis "L'analyse statistique des trajectoires" (see references).

It provides the datasets used in the examples in the handbook, as well as functions for :

* describing episodes in individual sequences (at least one episode, number of episodes, position of the start of the first episode)
* measuring association between domains in multidimensional sequence analysis
* heat maps of sequence data
* Globally Interdependent Multidimensional Sequence Analysis (GIMSA)
* smoothing sequences for index plots
* coding sequences for Qualitative Harmonic Analysis
* measuring stress from MDS factors
* symmetrical PLS


## Documentation

Please visit [https://nicolas-robette.frama.io/seqhandbook/](https://nicolas-robette.frama.io/seqhandbook/) for documentation


## installation

Execute the following code within `R`:

``` r
if (!require(devtools)){
    install.packages('devtools')
    library(devtools)
}
install_git("https://framagit.org/nicolas-robette/seqhandbook")
```

## References

Robette, Nicolas. *L'analyse statistique des trajectoires : Typologies de séquences et autres approches*. Nouvelle édition [en ligne]. Paris : Ined Éditions, 2021. Disponible sur Internet : <https://books.openedition.org/ined/16670>. ISBN : 9782733290507. DOI : https://doi.org/10.4000/books.ined.16670.
