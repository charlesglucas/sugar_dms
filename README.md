SUGAR D-MS
===

## Description
SUGAR DMS is an image recontruction and contour detection matlab package. It performs the minimization of the Discrete Mumford-Shah (D-MS) functional, which enforces constraints related to smoothing over image and sparsity over contours to obtain a piecewise smooth reconstructed image $u$ and sparse estimated contours $e$ from observed noisy images $z$:
$ \minimize_{u,e} $ $\Vert u - z \Vert_2$ + $ \beta \Vert (1-e) \odot Du \Vert_2 + \lambda \Vert e \Vert_2$

where $\odot$ denotes the component-wise product and $\beta$ > 0 and $\lambda > 0 are regularization parameters.

A Stein-like strategy providing optimal hyperparameters is designed, based on the minimization of an unbiased quadratic risk estimate. Efficient and automated minimization of the risk estimate relies on an unbiased estimate of the risk's gradient with respect to hyperparameters.

![alt text](http://perso.ens-lyon.fr/charles.lucas/images/DMSdenoising.svg)

## Requirements
The following matlab package is required: [GRANS0](https://gitlab.com/timmitchell/GRANSO/).

## References
  - [Lucas et al., 2021](https://arxiv.org/pdf/2109.13651.pdf)
  - [Foare et al., 2019](https://hal.archives-ouvertes.fr/hal-01782346/document)
  - [Deledalle et al., 2014](https://arxiv.org/pdf/1405.1164)
  
## Quick start
The basic syntax to run SUGAR D-MS is as follows:

```
% return optimal hyperparameters of D-MS
[Lambda,~] = bfgs_sugar_dms(image);

% return D-MS reconstructed image u and contour estimate e
[u,e,~] = DMS_2D(image,Lambda(1),Lambda(2));
```


The main parameters to take into account are:

  - `R`, the number of realizations of the Monte Carlo vector;
  - `sigma`, the noise level which is estimated by default.
    
Here is an example with non-default parameters:
```
param.R = 5; param.sigma = 0.1;
[Lambda,~] = bfgs_sugar_dms(image,param);
[u,e,~] = DMS_2D(image,Lambda(1),Lambda(2));
```
An example with simulated images can also be found in the `example` folder.
