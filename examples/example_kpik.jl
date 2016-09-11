######################################
# Example of kpik function
# Based on example_kpik.m on V. Simoncini's Website
######################################
workspace()
include("../src/matrixEqs.jl")
using matrixEqs

nh = 15
n = nh^2
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
A = sparse(-(kron(T,I)+kron(I,T)))
#A[7,1] = 1
srand(1234)
B = randn(n,2)
E = sparse(diagm(rand(n),0))
#LE=full(cholfact(E)[:L])
m=100
tol=1e-9
tolY=1e-12
Z,er2=kpik(A,B,E,tolY=tolY)

#using ProfileView
#kpik(A,B,E,tolY=tolY)
#Profile.clear()
#@profile kpik(A,B,E,tolY=tolY)
#ProfileView.view()


#norm(A*Z*Z'*E+E*Z*Z'*A'+B*B') > 1 && println(i)
println(norm(A*Z*Z'*E+E*Z*Z'*A'+B*B'))
