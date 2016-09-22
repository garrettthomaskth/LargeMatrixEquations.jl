# LME_Pack (Large Matrix Equations Package)
LME_Pack is a package that implements popular Numerical Methods for Large Scale Matrix Equations in Julia.

####Included methods for solving the following matrix equations
###Generalized Lyapunov Equation:
######Extended Krylov Method with Galerkin condition (kpik)
Based on kpik.m from V. Simoncini's [Website](http://www.dm.unibo.it/~simoncin/software.html). 
Resources:
V. Simoncini, 
A new iterative method for solving large-scale Lyapunov matrix equations 
SIAM J. Scient. Computing, v.29, n.3 (2007), pp. 1268-1288. 

######Adaptive Rational Krylov Subspaces Method for Lyapunov Equations
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
