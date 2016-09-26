# LME_Pack (Large Matrix Equations Package)
**LME_Pack** is a package that implements popular Numerical Methods for Large Scale Matrix Equations in Julia.

####Included methods for solving the following matrix equations
###Generalized Lyapunov Equation:
######Extended Krylov Method with Galerkin condition (kpik)
Based on kpik.m from V. Simoncini's [Website](http://www.dm.unibo.it/~simoncin/software.html). 
Resources:
V. Simoncini, 
A new iterative method for solving large-scale Lyapunov matrix equations 
SIAM J. Scient. Computing, v.29, n.3 (2007), pp. 1268-1288. 

######Adaptive Rational Krylov Subspaces Method (rksm) for Lyapunov Equations
Based on rksm.m from V. Simoncini's [Website](http://www.dm.unibo.it/~simoncin/software.html).
Resources:
V. Druskin and V. Simoncini, 
Adaptive rational Krylov subspaces for large-scale dynamical systems 
Systems & Control Letters, 60 (2011), pp. 546-560. 
Convexhull for Julia [GitHub](https://github.com/intdxdt/convexhull.jl).

######Low Rank Cholesky Factor Alternating Direction Implicit (LRCF-ADI) Method
Based on lp_lradi.m from LAYPACK 1.0, [Website](https://www.tu-chemnitz.de/sfb393/lyapack/).
Resources:
J.Li, F.Wang, and J.White.
An efficient Lyapunov equation based approach for generating
reduced-order models of interconnect.
Proceedings of the 36th IEEE/ACM Design Automation Conference,
New Orleans, LA, 1999.

#Example
```julia
using LME_Pack

n = 10
A = -diagm(ones(nh)*2)+diagm(ones(nh-1),-1)+diagm(ones(nh-1),1)
B = [1:n (n+1):(2*n)]

Zkpik,er2=kpik(A,B)
Zrksm,resnorm=rksm(A,B,E)
Zadi,flag,res=lp_lradi(A,B)

# Example of ploting the backwards error from kpik
using PyPlot
title("Plot of Scaled Residual")
ylabel("Residual")
xlabel("Iteration")
semilogy(er2, color="red", linewidth=2.0, linestyle="--")
show()
```

```julia


 ```
