# <div align="center">SUGAR D-MS</div> 

## Reference

> **[Charles-Gérard Lucas](https://perso.ens-lyon.fr/charles.lucas), [Barbara Pascal](https://bpascal-fr.github.io), [Nelly Pustelnik](http://perso.ens-lyon.fr/nelly.pustelnik/), [Patrice Abry](https://perso.ens-lyon.fr/patrice.abry),**
*Hyperparameter selection for the Discrete Mumford-Shah functional,* 
Preprint. [Download](https://arxiv.org/pdf/2109.13651.pdf)

## Description
SUGAR DMS is an image recontruction and contour detection matlab package. It performs the minimization of the Discrete Mumford-Shah (D-MS) functional, which enforces constraints related to smoothing over image and sparsity over contours to obtain a piecewise smooth reconstructed image $u$ and sparse estimated contours $e$ from observed noisy images $z$:
$$\min_{u,e} \Vert u - z \Vert_2 +  \beta \Vert (1-e) \odot Du \Vert_2 + \lambda \Vert e \Vert_1$$

where $\odot$ denotes the component-wise product, D is a discrete difference operator, and $\beta > 0$ and $\lambda > 0$ are regularization parameters.

A Stein-like strategy providing optimal hyperparameters $\beta$ and $\lambda$ is designed, based on the minimization of an unbiased quadratic risk estimate. Efficient and automated minimization of the risk estimate relies on an unbiased estimate of the risk's gradient with respect to hyperparameters.

<p align="center">
  <img width="300" src="http://perso.ens-lyon.fr/charles.lucas/images/DMSdenoisingIllustration.svg">
</p>

Here is a comparison with the state-of-the-art SUGAR T-ROF (see [Cai and Steidl, 2013](https://page.math.tu-berlin.de/~steidl/PDFs/CaiSte13.pdf) for T-ROF and [Deledalle et al., 2014](https://arxiv.org/pdf/1405.1164.pdf) for SUGAR) run using [gsugar](https://github.com/bpascal-fr/gsugar) on a [BSD69 dataset](https://paperswithcode.com/dataset/bsd) image:
<p align="center">
  <img width="1200" src="http://perso.ens-lyon.fr/charles.lucas/images/SUGARDMSresults.png">
</p>

## Recommendation
This toolbox is designed to work with [**Matlab 2020b**](https://fr.mathworks.com/products/new_products/release2020b.html).

## Requirement
The following matlab package is required: [GRANS0](https://gitlab.com/timmitchell/GRANSO/).

## Related works
  - [Foare et al., 2019](https://hal.archives-ouvertes.fr/hal-01782346/document)
  - [Deledalle et al., 2014](https://arxiv.org/pdf/1405.1164)
  
## Quick start
The basic syntax to run SUGAR D-MS is as follows:

```
[Lambda,~] = bfgs_sugar_dms(image); 
[u,e,~] = DMS_2D(image,Lambda);
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
