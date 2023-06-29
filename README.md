# svylme
Mixed models for complex surveys

This package fits linear mixed models to data from complex surveys, by maximising a weighted pairwise likelihood

## Advantages

It works (gives consistent estimates of the regression coefficients and variance components) for **any** linear mixed model and **any** design, without any restrictions on the sampling units
and model clusters being related. For example, you could sample on home address but fit a model clustering on school.

The implementation allows for correlated random effects such as you get in quantiative genetics

## Disadvantages

Some loss of efficiency compared to just fitting a design-based linear model (if you don't care about the variance components)

There isn't (yet) an analog of the BLUPs of random effects

If your sampling units and model clusters are the same, and your design isn't too strongly informative, you can likely get more precise estimates of the variance components with
sequential pseudolikelihood as implemented in Stata or Mplus. 
