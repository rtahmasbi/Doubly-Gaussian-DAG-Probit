# Bayesian Causal Inference in Doubly Gaussian DAG-probit Models

This repo is related to my paper: https://arxiv.org/abs/2304.05976

`R` codes are under the `R` directory.


## summary
We consider modeling a binary response variable together with a set of covariates for two groups under observational data. The grouping variable can be the confounding variable (the common cause of treatment and outcome), gender, case/control, ethnicity, etc. Given the covariates and a binary latent variable, the goal is to construct two directed acyclic graphs (DAGs), while sharing some common parameters. The set of nodes, which represent the variables, are the same for both groups but the directed edges between nodes, which represent the causal relationships between the variables, can be potentially different. For each group, we also estimate the effect size for each node. We assume that each group follows a Gaussian distribution under its DAG. Given the parent nodes, the joint distribution of DAG is conditionally independent due to the Markov property of DAGs. We introduce the concept of Gaussian DAG-probit model under two groups and hence doubly Gaussian DAG-probit model. To estimate the skeleton of the DAGs and the model parameters, we took samples from the posterior distribution of doubly Gaussian DAG-probit model via MCMC method. We validated the proposed method using a comprehensive simulation experiment and applied it on two real datasets. Furthermore, we validated the results of the real data analysis using well-known experimental studies to show the value of the proposed grouping variable in the causality domain.


# Causality
To learn more about causality, check my repo at https://github.com/rtahmasbi/Causality
