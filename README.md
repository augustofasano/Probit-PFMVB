Asymptotically Exact Variational Bayes for High-Dimensional Probit Models
=========================================================================

This repository is associated with the article (ADD LINK) Asymptotically
Exact Variational Bayes for High-Dimensional Probit Models. The **key
contribution of this paper is outlined below**.

*In this article we propose a novel variational approximation for the
posterior distribution of the coefficients in a probit regression with
Gaussian priors. Our method leverages a representation with global and
local variables but, unlike classical mean-field assumptions, we
crucially avoid a fully factorized approximation, and instead rely on a
variational family in which only the joint density of the local
variables is factorized.*

This repository provides **codes and tutorials to implement the
inference methods associated with such a new result**. In particular,
the focus is on the case-study outlined in Section 3 of the paper. The
complete tutorial can be found in the file
<tt>AlzheimerTutorial.html</tt>.

This tutorial focuses on a large *p* and small-to-moderate *n*
biomedical application. The dataset is obtained by a modification of an
original dataset arising from a study of the Washington University to
determine if biological measurements from cerebrospinal fluid are useful
in modeling and predicting early stages of Alzheimer’s disease
([Craig-Shapiro et al.,
2011](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0018850)).
Full details on the preprocessing are available in the tutorial file.

The goal is to compare the **Algorithm 2** proposed in the paper—which
provides a partially-factorized mean-field (PFM) approximation the
unified skew-normal posterior—with state-of-the-art mean-field (MF)
approximations in [Consonni and Marin
(2007)](https://www.sciencedirect.com/science/article/pii/S0167947306003951),
reported in **Algorithm 1**, using i.i.d. samples from the exact unified
skew-normal posterior distribution in [Durante
(2019)](https://arxiv.org/abs/1802.09565) as a benchmark for the
accuracy.

**The methods involved in the analysis are implemented though
functions**, that can be found in the file
<tt>functionsVariational.R</tt>. A **tutorial explaining in details the
usage of such functions** is also made available, for additional
clarity, in the file <tt>functionsTutorial.html</tt>.

All the analyses are performed with a MacBook Pro (OS Mojave, version
10.14.6, Processor 2.7 GHz Intel Core i5, RAM 8 GB), using an <tt>R</tt>
version 3.6.1.

IMPORTANT: Although a seed is set at the beginning of each sampling
scheme, the final output reported in Tables and Figures of
AlzheimerTutorial.html may be subject to slight variations depending on
which version of the <tt>R</tt> packages has been used in the
implementation of the code. This is due to possible internal changes of
certain functions when the package version has been updated. **However,
the magnitude of these minor variations is negligible and does not
affect the final conclusions**.
