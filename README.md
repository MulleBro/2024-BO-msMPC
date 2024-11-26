# Automatic design of robust model predictive control of a bioreactor via Bayesian Optimization

Accompanying repository for our conference paper "Automatic design of robust model predictive control of a bioreactor via Bayesian Optimization".
For questions about these materials we kindly ask you to send an e-mail at tobias.brockhoff@tu-dortmund.de.

## Abstract

Model predictive control (MPC) is an advanced control strategy that can deal
with general nonlinear systems and constraints but relies on accurate predictions given by
a dynamic model. To satisfy constraints and improve performance despite imperfect models,
robust MPC methods can be formulated. Multi-stage MPC is a robust MPC method based
on the formulation of scenario trees. However, the resulting optimization problems can be
very large, as the number of scenarios considered in the tree results from the combinations
of all possible uncertainties. For systems when the number of uncertainties is large, as it is the
case in bioprocesses, the optimization problems become rapidly intractable. To solve this issue,
heuristics are typically used to select the most relevant uncertain parameters and their range
of uncertainty. In this paper, we propose a two-step approach to obtain a systematic design of
multi-stage MPC controllers: First, the key uncertain parameters are extracted based on the
parametric sensitivities. Second, Bayesian Optimization is employed for tuning of the range of
uncertainties. The approach is applied to a bioreactor simulation study. The proposed approach
can avoid constraint violations that are otherwise obtained by standard MPC while being less
conservative than a manually-tuned robust controller.

## About this repository
All the code to generate the results are provided as well as the generated results themselves.
The documentation is currently under construction
