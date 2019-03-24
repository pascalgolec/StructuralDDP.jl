# DiscreteDynamicProgramming

This package solves and simulates discrete dynamic choice models with value function iteration.

<img src="https://latex.codecogs.com/svg.latex?\alpha^2" title="\alpha^2" />

## Problem definition

Each discrete dynamic problem has can be written in the fol

Each discrete dynamic problem has

- a flow reward as a function of the state variables and choice(s)
- transition equation(s) for state variables
	- --> for calculating probability transition matrix
	- an also supply own prob transition matrix
- parameters
- a tuple of vectors of state variables
	- tStateVectors
	- the choice state variables need to come first
- a tuple of vectors of choice variables
	- tChoiceVectors

## Solver assumptions

- state variables are always forced to stay within bounds
- for `intermediate` and `separable`, can choose directly the state variable
- order of state variables is first the endogenous state variable, then exogenous

## Reserved parameters

- discount factor beta
- solver options
