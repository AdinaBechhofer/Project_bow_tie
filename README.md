# Project_bow_tie
## Course project for 6.336
### Links to writeups:
https://www.overleaf.com/9353475874xbfcgntrjpzb
< insert link to project reports>

### Code Documentation

#### eval_f3
This function evaluates ![equation](https://latex.codecogs.com/gif.latex?\frac{dx}{dt}&space;=&space;f(x,&space;p,&space;u,&space;b))
Usage: eval_f3 is called `eval_f3(x,p, u, b, t)` where the t is an optional argument. Likewise, it can be called `eval_f3(x, p, u, b)`

x: Vector of length = number of state variables.

p: Struct with the relevant parameters. Must contain:
    
    * p.Area
    
    * p.Beta
    
    * p.Distance
    
    * p.workFunction
    
    * p.Radius
    
    * p.taby
    
    * p.invC
    
    * p.CG

u: Vector of inputs 

b: Matrix of length = length(x) and width = length(u)

#### fjbowtie
This function returns the f and J for eval_f3. 

Usage: `fjbowtie(x,p,u,b,t)` or `fjbowtie(x,p,u,b)`. The t is optional 

u: Vector. Must be precomputed. 

#### FiniteDifferenceJacobian_t
This function computes the finite difference Jacobian of the handle function given to it

Usage: `FiniteDifferenceJacobian_t(f,x,p,u,b,t)` or `FiniteDifferenceJacobian_t(f,x,p,u,b)`. The `t` is optional. 

#### ForwardEuler_t
This function computes the forward Euler simulation. 

Usage: `X = ForwardEuler_t(fhand, x0,p,u,b,t)` or `[X, t] =  ForwardEuler_t(fhand, x0,p,u,b, t_start, t_stop, del_t)` or `[X, t] =  ForwardEuler_t(fhand, x0,p,u,b, 'dynamic', t_start, t_stop)`

