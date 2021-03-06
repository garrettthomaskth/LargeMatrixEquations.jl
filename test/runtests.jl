workspace()
include("../src/LargeMatrixEquations.jl")
using LargeMatrixEquations
using Base.Test

nh = 10
n = nh^2#nh^2
T = diagm(ones(nh)*2)-diagm(ones(nh-1),-1)-diagm(ones(nh-1),1)
I = eye(nh,nh)
A = -(diagm(ones(n)*2)-diagm(ones(n-1),-1)-diagm(ones(n-1),1))#-(kron(T,I)+kron(I,T))
srand(1234)
B = ones(n,1)

E = diagm(ones(n),0)
E[1,1]=2
E[n,n]=3

Z,er2=try
        kpik(A,B,E,infoV=false)
      catch
        error("kpik failed with full matrices")
      end
@test norm(A*Z*Z'*E + E*Z*Z'*A' + B*B') < 0.001


Z,er2=try
        kpik(sparse(A),B,sparse(E))
      catch
        error("kpik failed with sparse matrices")
      end

@test norm(A*Z*Z'*E + E*Z*Z'*A' + B*B') < 0.001


Z,resnorm=try
              rksm(A,B,E,infoV=false)
          catch
              error("rksm failed with full matrices")
          end
@test norm(A*Z*Z'*E + E*Z*Z'*A' + B*B') < 0.001

Z2,resnorm2=try
            rksm(sparse(A),B,sparse(E))
          catch
              error("rksm failed with sparse matrices")
          end
@test norm(A*Z2*Z2'*E + E*Z2*Z2'*A' + B*B') < 0.001


Z,flag,res =try
              lp_lradi(A,B)
            catch
              error("lp_lradi failed with full matricies")
            end
@test norm(A*Z*Z' + Z*Z'*A' + B*B') < 0.001

Z,flag,res =try
              lp_lradi(sparse(A),B,infoV=false)
            catch
              error("lp_lradi failed with sparse matricies")
            end
@test norm(A*Z*Z' + Z*Z'*A' + B*B') < 0.001
