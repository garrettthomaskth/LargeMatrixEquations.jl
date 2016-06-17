cd("/Users/garrettthomas/matrixEqs")
include("kpikSimp.jl")
nh = 3
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
B = randn(n,1)
#B = Float64(n,1)
#for i = 1:n
#  B[i,1]= (i+1)/10
#end
m=100
tol=1e-9
tolY=1e-12
kpikSimp(A,B,m,tol,tolY)
