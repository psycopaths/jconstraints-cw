# jConstraints-cw #
Implementation of the Concolic Walk algorithm in jConstraints. The algorithm is described in [Solving complex path conditions through heuristic search on induced polytopes](http://dl.acm.org/citation.cfm?doid=2635868.2635889) by Peter Dinges and Gul Algha.

## Installation ##
First, install [jConstraints](https://bitbucket.org/teamcoco/jconstraints).
Then, ```mvn install``` in the root of the jconstraints-cw folder.

To enable the solver, put in your jpf file
```
#!txt
symbolic.dp=cwsolver
```

## Configuration ##
Below is the jpf configurations for the Concolic Walk solver and their default values:

```
#!txt
concolic_walk.linear_solver=z3
concolic_walk.iterations_per_constraint=150
concolic_walk.neighbors_generated_per_iteration=10
concolic_walk.tabu_iterations_per_variable=0.5f
concolic_walk.min_tabu_iterations=3
concolic_walk.enable_seeding=true
concolic_walk.enable_bisection=true
concolic_walk.randomization_radius_exponent=12
```

Explanation of the options:

### linear_solver ###
Which solver to use for solving the linear constraints. It can be set to the string value of *any* installed jConstraints-solver.

### iterations_per_constraint ###
Number of algorithm iterations added for each constraint in the path condition. The paper speaks of *steps* instead of iterations and denotes this constant as **I**. 

### iterations_generated_per_iteration ###
Number of random neighbor RealVector points to generate per iteration. The paper calls points *environments* and denotes this constant as **R**.  

### tabu_iterations_per_variable ###
Multiplier for computing how many iterations a variable is tabu if changing it failed to yield a better point. The paper uses this constant and min_tabu_iterations to compute the number of tabu iterations **T**.

### min_tabu_iterations ###
Minimum number of iterations that a variable is tabu if changing it failed to yield a better point. See tabu_iterations_per_variable  above.

### enable_seeding ###
Enable trying the hard-coded SEEDED_VALUES for the changed variable?

### enable_bisection ###
Enable estimating additional neighbors that *should* satisfy one of the constraints violated by the current variable if the constraint was linear?  See the discussion of Finding Neighbors in section 3.3 of the paper.


## Limitations ##
See some of the test cases in the test suites. In summary:

* Having a constraint composed of **only** nonlinear conjuncts, cannot be solved by CW, since we cannot build a polytype. A solution could be to try random inputs if a polytype cannot be constructed. In fact, if it wasn't for floating point imprecision, we could, by randomly assigning values to the variables solve a constraint such as ```('x' >= 1.2) && ('tan'('x') == ('sin'('x') / 'cos'('x'))))``` (the identity)

* If non-linear constraints are present (and therefore actually applying the concolic walk), the solver will **not** be able to conclude UNSAT; the result can only be DONT_KNOW or SAT. This should maybe be investigated further since it may be due to adaption of the original algorithm

* That example demonstrates another limitation: because the nonlinear constraints (and native methods) **are** evaluated, imprecisions in floating point computations can result in false conclusions e.g. even though we can come up with a correct assignment for the identity  ```('x' >= 1.2) && ('tan'('x') == ('sin'('x') / 'cos'('x'))))```, CW will likely give DONT_KNOW.

## Future Work ##
* Combine with other solvers supporting nonlinear constraints e.g. coral. For example, if the constraint is only composed of nonlinear conjuncts (which CW cannot solve -- see above), resort to just use coral.

## Licensing ##
The implementation of the Concolic Walk algorithm is available under the MIT license. Copyright (c) 2014 the University of Illinois Board of Trustees. All Rights Reserved.

The Concolic Walk algorithm is available on [the GitHub page](https://github.com/osl/concolic-walk).