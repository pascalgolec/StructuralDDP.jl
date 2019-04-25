# StructuralDDP

This package defines, solves and simulates discrete dynamic optimization problems with value function iteration fast. It features different options to add structure to the problem and provide information about its properties. A wide range of problems can be solved very fast this way.

The following figure compares the speed of standard solution methods (value function iteration using a state-action representation) with the structured approach for solving a capital investment model with one or two types of capital.

![alt text](benchmark/compare.svg "Benchmarking acceleration")

In the structured approach we specified that only some of the state variables should be integrated when calculating expectations, that the optimal choices are monotone in some of the states and that the value function is concave in some of the states and finally that part of the reward matrix should be pre-built before calling the solver.

The package also allows problems to be defined using their mathematical formulation with minimal additional input. There is no need to calculate any matrices for example. One only has to provide a reward function, a transition function, as well as state and action spaces.
