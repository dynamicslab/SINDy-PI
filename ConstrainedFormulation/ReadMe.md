# Instruction

This folder constains the example of constrained formulation of SINDy-PI. We are interested in solving the problem

$$
\min_{\boldsymbol{\Xi}} \quad  \|\boldsymbol{\Theta}(\boldsymbol{X},\boldsymbol{\dot{X}})-\boldsymbol{\Theta}(\boldsymbol{X},\boldsymbol{\dot{X}})\boldsymbol{\Xi}\|_2+\beta\|\boldsymbol{\Xi}\|_0
\\
\quad \text{s.t. }\quad \text{diag}(\boldsymbol{\Xi})=\boldsymbol{0}.
$$

Please open the Main.m file to run this example.
![](Images/SINDy_PI_Constrained.jpg)
# Prerequisite

In order to run this file, you'll need to instal the [CVX]([http://cvxr.com/cvx/](http://cvxr.com/cvx/)) optimization package. Please click the link for details.
