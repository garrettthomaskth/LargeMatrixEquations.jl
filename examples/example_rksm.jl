######################################
# Example of rksm function
# Based on example_rksm.m on V. Simoncini's Website
######################################
workspace()
include("../src/matrixEqs.jl")
using matrixEqs
nh=30;
T = zeros(nh,nh)
T[1,1] = 2
T[1,2] = -1
T[nh,nh-1] = -1
T[nh,nh] = 2
for i in 2:(nh-1)
  T[i,i-1] = -1
  T[i,i] = 2
  T[i,i+1] = -1
end
T[2,1] = 4
I=eye(nh);
A=-(kron(T,I)+kron(I,T));
n=nh^2;
srand(123)
E = diagm(rand(n),0)
EL=cholfact(E)[:L]
B=randn(n,2);
m=100;
tol=1e-9;
tolY=1e-12;
ch=true;
Z,resnorm=rksm(A,E,EL,B)#,m,tol,s1,emax,ch,tolY);
print(norm(A*Z*Z'*E+E*Z*Z'*A'+B*B'))
