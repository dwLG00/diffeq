# Differential Equations Modeling

Mainly functions/programs for ease of modeling differential equations and their solutions

## Artifact 1

In `approximations.py`, there are two ways of "estimating" the solution of a differential equation, given a starting point.

The `solve` function and its wrappers utilize an Euler approximation to estimate the solution. Given a small dx, it iteratively calculates the "next" point and its derivative at that point. The `picard` function instead uses Picard iteration to compute an approximation for the solution.

The difference in the two approaches comes down to soundness. In the euler approximation, in the case where there exists more than one solution past some threshold, the approximation will derive one of the solutions but not the other. In the Picard iteration approach, instead, we know that the solution found must be unique for some interval.
