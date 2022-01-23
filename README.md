# fun2optimize

Fun2Optimize is a project realized in MATLAB which helps you understand tunning parameters of various hand-coded optimization methods. 

Explore various optimization methods and play around with tunning parameter values to see their influence.

Main.mlx shows implementation on optimization problem with nonlinear equality constraint - rectangle to be inscribed in the ellipse - that can be reduced to an unconstrained problem using the elimination method.

Currently supported optimization methods:
- Gradient Method
- Gradient Method with Backtracking
- Newton Method
- Luus-Jaakola Method
- Simulated Annealing

## Getting Started
1. Formulate the optimization problem as function handle
	```
	f = @(x) -4*x*sqrt(1-x^2/4)
	```

2. Find the solution of the optimization problem

## TODO

- [ ] Add Python implementation
- [ ] Add multiparameter optimization support
