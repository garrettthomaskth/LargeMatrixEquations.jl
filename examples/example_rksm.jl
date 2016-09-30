######################################
# Example of rksm function
# Based on example_rksm.m on V. Simoncini's Website
######################################
workspace()
include("../src/LME_Pack.jl")
using LME_Pack

nh = 20
n = nh^2
T = diagm(ones(nh)*2)-diagm(ones(nh-1),-1)-diagm(ones(nh-1),1)
I = eye(nh,nh)
A = sparse(-(kron(T,I)+kron(I,T)))

srand(123)

E = eye(n)#diagm(rand(n),0)

B=randn(n,2)

m=100
tol=1e-9
tolY=1e-12
ch=true
Z,resnorm=rksm(A,B)
print(norm(A*Z*Z'*E+E*Z*Z'*A'+B*B'))
