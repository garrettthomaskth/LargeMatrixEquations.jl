# LargeMatrixEquations
**LargeMatrixEquations** is a package that implements popular Numerical Methods for Large Scale Matrix Equations in Julia.

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

###Download
git clone https://github.com/garrettthomaskth/LargeMatrixEquations.git

###Example
```julia
using LargeMatrixEquations

n = 1000
A = -diagm(ones(n)*2)+diagm(ones(n-1),-1)+diagm(ones(n-1),1)
B = ones(n,1)

Zkpik,er2=kpik(A,B,tol=1e-14)
println(norm(A*Zkpik*Zkpik' + Zkpik*Zkpik'*A' + B*B'))

Zrksm,resnorm=rksm(A,B,tol=1e-14)
println(norm(A*Zrksm*Zrksm' + Zrksm*Zrksm'*A' + B*B'))


Zadi,flag,res=lp_lradi(A,B)
println(norm(A*Zadi*Zadi' + Zadi*Zadi'*A' + B*B'))

println("LYAP")
C=B*B'
tic()
X=lyap(A,C*1.)
println(toq())
println(norm(A*X + X*A' + B*B'))
```
This example shows how one can compare the preformence of the three methods included in LargeMatrixEquations with the built in Lyapunov solver, ```lyap```.

```julia
using PyPlot
title("Plot of KPIK Backwards Error")
ylabel("Error")
xlabel("Iteration")
semilogy(er2, color="red", linewidth=2.0, linestyle="--")
show()
```

The second part of this example shows how one can use the package PyPlot to visualize the information returned from ```kpik```. This idea is of course also applicable to the other functions. 
