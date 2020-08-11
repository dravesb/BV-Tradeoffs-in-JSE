# Bias-Varince-Tradeoffs-in-Joint-Spectral-Embeddings

Latent position models and their corresponding estimation procedures offer a statistically principled
paradigm for multiple network inference by translating multiple network analysis problems to familiar
tasks in multivariate statistics. Latent position estimation is a fundamental task in this framework yet
most work focus only on unbiased estimation procedures. We consider the ramifications of utilizing biased
latent position estimates in subsequent statistical analysis in exchange for sizable variance reductions in
finite networks. We establish an explicit bias-variance tradeoff for latent position estimates produced by
the omnibus embedding of Levin et al. (2017) in the presence of heterogeneous network data. We reveal
an analytic bias expression, derive a uniform concentration bound on the residual term, and prove a
central limit theorem characterizing the distributional properties of these estimates. These explicit bias
and variance expressions enable us to show that the omnibus embedding estimates are often preferable to
comparable estimators with respect to mean square error, state sufficient conditions for exact recovery in
community detection tasks, and develop a test statistic to determine whether two graphs share the same
set of latent positions. These results are demonstrated in several experimental settings where community
detection algorithms and hypothesis testing procedures utilizing the biased latent position estimates are
competitive, and oftentimes preferable, to unbiased latent position estimates.

Check out the paper pre-print here: <https://arxiv.org/abs/2005.02511>

[test](one_dim_BV_tradeoff/figures/1d_mse.pdf)
