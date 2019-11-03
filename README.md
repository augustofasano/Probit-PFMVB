# Asymptotically Exact Variational Bayes for High-Dimensional Binary Regression Models

This repository is associated with the article [Fasano, Durante and Zanella (2019). *Asymptotically Exact Variational Bayes for High-Dimensional Binary Regression Models*]. The **key contribution of this paper is outlined below**.

> In this article we propose a novel variational approximation for the posterior distribution of the coefficients in a probit regression with Gaussian priors. Our method leverages a representation with global and local variables but, unlike classical mean-field assumptions, we crucially avoid a fully factorized approximation, and instead rely on a variational family in which only the joint density of the local variables is factorized.

This repository provides **codes and tutorials to implement the inference methods associated with such a new result**. In particular, the focus is on the large *p* application to Alzheimer’s data outlined in Section 3 of the paper (see also [Craig-Shapiro et al., 2011](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0018850)). The complete tutorial can be found in the file [`AlzheimerTutorial.md`](https://github.com/augustofasano/Probit-PFMVB/blob/master/AlzheimerTutorial.md) where we also provide details to pre-process the original dataset available in the `R` package `AppliedPredictiveModeling`. 

The goal of this analysis is to compare the performance of the proposed **partially-factorized mean-field (PFM) approximation** relative to those state-of-the-art competitors which were feasible in this application. These include the classical **mean-field (MF) variational approximation** (see [Consonni and Marin, 2007](https://www.sciencedirect.com/science/article/pii/S0167947306003951)) and **Monte Carlo inference based on i.i.d. from the exact unified skew-normal posterior distribution** derived by [Durante (2019)](https://doi.org/10.1093/biomet/asz034). The latter serves also a benchmark to study the accuracy of the approximate methods. See Section 2 in the article for details. We also tried to implement **Hamiltonian Monte Carlo methods** (`R` package `rstan`) and **expectation-propagation** (`R` package `EPGLM`), but these algorithms were impractical. Hence, we will not focus on such schemes in this repositiory.

**The functions to implement the above methods** can be found in the `R` source file [`functionsVariational.R`](https://github.com/augustofasano/Probit-PFMVB/blob/master/functionsVariational.R), and a **tutorial explaining in detail the usage of these functions** is available in the file [`functionsTutorial.md`](https://github.com/augustofasano/Probit-PFMVB/blob/master/functionsTutorial.md).

All the analyses are performed with a **MacBook Pro (OS Mojave, version 10.14.6, Processor 2.7 GHz Intel Core i5, RAM 8 GB)**, using an `R` version **3.6.1**.

**IMPORTANT**: Although a seed is set at the beginning of each sampling scheme, the final output reported in Tables and Figures of [`AlzheimerTutorial.md`](https://github.com/augustofasano/Probit-PFMVB/blob/master/AlzheimerTutorial.md) may be subject to slight variations depending on which version of the `R` packages has been used in the implementation of the code. This is due to possible internal changes of certain functions when the package version has been updated. **However, the magnitude of these minor variations is negligible and does not affect the final conclusions**.
