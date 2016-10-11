######################################
# Example of kpik function
# Based on example_kpik.m on V. Simoncini's Website
######################################
workspace()
include("../src/LargeMatrixEquations.jl")
using LargeMatrixEquations

nh = 20
n = nh^2
T = diagm(ones(nh)*2)-diagm(ones(nh-1),-1)-diagm(ones(nh-1),1)
I = eye(nh,nh)
A = sparse(-(kron(T,I)+kron(I,T)))

srand(1234)
B = randn(n,2)
E = eye(n)#sparse(diagm(rand(n),0))

m=100
tol=1e-9
tolY=1e-12
Z,er2=kpik(A,B,tolY=tolY)

println(norm(A*Z*Z'*E+E*Z*Z'*A'+B*B'))
