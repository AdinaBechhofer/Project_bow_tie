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

Usage: `X = ForwardEuler_t(fhand, x0,p,U,b,t)` or `[X, t] =  ForwardEuler_t(fhand, x0,p,U,b, t_start, t_stop, del_t)` or `[X, t] =  ForwardEuler_t(fhand, x0,p,U_func,b, 'dynamic', t_start, t_stop)`. Where `U` is the time dependant input matrix of size `length(x0)*length(tvec)`

Other Usages: To perform forward Euler on a linearized system with coefficients [A, B], `X = ForwardEuler_t(fhand,x0, p, U, b, 'linearmodel', t, A, B)`
 

#### FJFTrap 
This function computes the f and J for the trapezoidal method 

Usage: `FJFTrap(x, p, u, b, t, gamma, dt, fJhand)`


#### newtonNd
This function does the newton method for an N dimensional vector. 

Usage: `newtonNd(fJfhand,x0,p, u, b, t)` or `newtonNd(fJfhand,x0,p, u, b)` for a normal non-linear function computation. 

For trapezoidal method solving, use `newtonNd(TrapHand,x0,p, u, b, t, gamma, dt, integrand)` where `TrapHand` is `FJFTrap` and `integrand` is the fJ function handle for the functin you're interested in integrating. 


 #### TrapMethod 
 This function performs the trapezoidal time integration method for an input function. 

 Usage: `TrapMethod(x0,p,u,b, fJfhand, tvec)` where `fJhand` is the function you want to integrate. 
 

### Other misc. code

#### linearization

This function returns the linearized coefficients of the system [A, B] for linearization about input state x0 AND bias point u0; or just linearization about the input state x0.

Usage: To linearize only about x0, `linearization(f,x0,p,u,b,t,'onlyx0')`. The `t` is neccessary, but `u` and `t` are just dummy variables.
To linearize about both x0 and u0, `linearization(f,x0,p,u,b,t)` or `linearization(f,x0,p,u,b)`. The `t` is optional.


