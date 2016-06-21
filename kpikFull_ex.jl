######################################
# Script to test kpikFull.jl
# Based on example_kpik.m in the Davide folder
######################################

cd("/Users/garrettthomas/matrixEqs")
include("kpikFull.jl")
nh = 30
n = nh^2
# Create the Matrix T (find a more sophisticated way to do this)
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
I = eye(nh,nh)
A = -(kron(T,I)+kron(I,T))
B = randn(n,2)
E = diagm(rand(n),0)
LE=full(cholfact(E)[:L])
m=100
tol=1e-9
tolY=1e-12
kpikFull(A,E,LE,B,m,tol,tolY)
