[![Build Status](https://travis-ci.org/lruthotto/OptimTools.jl.svg?branch=master)](https://travis-ci.org/lruthotto/OptimTools.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/61rviuke750imsc1?svg=true)](https://ci.appveyor.com/project/lruthotto/optimtools-jl)
[![Coverage Status](https://coveralls.io/repos/lruthotto/OptimTools.jl/badge.svg)](https://coveralls.io/r/lruthotto/OptimTools.jl)




OptimTools.jl
=========================

&copy; 2015 [Lars Ruthotto](http://www.mathcs.emory.edu/~lruthot/). Released under the [MIT License](https://github.com/lruthotto/KrylovMethods.jl/blob/master/LICENSE).

Julia implementation of some tools for numerical optimization. ``OptimTools`` is primarily designed as a teaching tool. The first priority is to provide simple and easy-to-adapt implementations. The second priority is computational efficiency and testing.

## Installation

At the Julia prompt, type
```
julia> Pkg.clone("https://github.com/lruthotto/OptimTools.jl.git")
julia> Pkg.test("OptimTools")
```
## Quick Example
To minimize Rosenbrock's function using Newton's method, type

```
using OptimTools

f(x)   = (1-x[1])^2 + 100*(x[2]- x[1]^2).^2
df(x)  = [
            400*x[1]^3 - 400*x[1]*x[2] + 2*x[1]-2; 
            200*x[2]-200*x[1]^2  
         ]
d2f(x) = [
            1200*x[1]^2-400*x[2]+2      -400*x[1]; 
                         -400*x[1]           200
         ]

xnt, = newton(f,df,d2f,zeros(2))
```



