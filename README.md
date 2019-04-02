# DiscreteDynamicProgramming

This package solves and simulates discrete dynamic choice models with value function iteration fast. By taking advantage of properties of the problem, a wide array of problems can be solved very fast without requiring parallelization.

The following figure shows the time time to solve a model with one choice varaible (neoclassical) and two choice variables (intangibles) when using the acceleration methods in the toolbox versus not.

![alt text](benchmark/compare.svg "Benchmarking acceleration")

# To Do

- [ ] streamline model constructor. Put myself into shoes of user, try to create a model. End user shouldn't have to deal with types
- [ ] get rid of DDM, use DDP
- [ ] make tChoiceVectors optional (only if use SA)
    - in the problem constructor, make optional
    - remove it in separable and intermediate solver functions and other functions
    - change name of bEndogStatevars, to choicevars
- [ ] write documentation
- [ ] implement plotrecipe
- [ ] think about whether should include accountingvars() into this package
- [ ] :F is when the reward function has discontinuities
    - do later when design is fixed, because it may depend on it..
- [ ] firm exit
