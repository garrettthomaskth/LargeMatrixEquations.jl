# Julia code for K-PIK (Krylov-plus-inverted-Krylov)
# Based on kpik.m avalible from V. Simoncini's website
#
# Description from Simoncini:
#
# Approximately solve
#
#       A X E + E X A' + BB' = 0
#
# by means of the extended Krylov subspace method
# Input
#  A   coeff matrix, A < 0
#  E   coeff matrix, spd,
#  LE  lower triang factor of coeff matrix
#  B   factor of rhs,   nxk matrix with k << n
#  m   max space dimension, say sqrt(size(A))
#  tol stopping tolerance, with stopping criterion
#          ||LE\A X LE  + LE' X A'/LE'-LE\BB'/LE'||
#          ----------------------------------------  < tol
#      ||LE\BB'/LE'|| + ||E^{-1}|| ||A|| ||LE'X LE ||
#      computed in a cheap manner
#
#  Output:
#  Z   solution factor   X = Z Z'
#  er2 history of scaled residual, as above
#
#
# Comments:
# * The projected solution is computed at each iteration
#   As an alternative, a periodic computation could be considered.
# * This code performs a factorization of A. As an alternative,
#   iterative solves could be considered.
#
#  Please contact V. Simoncini for any problem you may encouter when
#  running the code
#
# If you use this code, please cite the following article:
#
# V. Simoncini
# A new iterative method for solving large-scale Lyapunov matrix equations,
# SIAM J.  Scient. Computing, v.29, n.3 (2007), pp. 1268-1288.
#
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
#FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
#COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
#IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

module kpikGL

export kpikFull

function kpikFull(A,E,LE,B,m=100,tol=1e-9,tolY=1e-12)
  @assert(isdefined(:vecnorm),"Your julia version is too old. vecnorm() not defined")
  tic();

  rhs=LE\B;
  nrmb=vecnorm(rhs)^2;
  nrma=vecnorm(A);
  sqrt2=sqrt(2);
  er2=zeros(m,1);

  n,sh=size(rhs);

  Y=[]
  odds=[]
  er2=[]

  if (vecnorm(E-speye(n))>1e-14)
    condestE=cond(E);
    singE=condestE/vecnorm(E);
  else
    singE=1;
  end

  if norm(A-A',1)<1e-14
     UA = chol(-A)
     LA = -UA'
     println("A sym. Completed Chol factorization\n")
     k_max =2
   else
     LA, UA=lu(A)
     println("A nonsym. Completed LU factorization\n")
     k_max = m
   end

   s=2*sh;
   rhs1=LE'*(UA\(LA\(LE*rhs)));

   # Orthogonalize [B,A^{-1}B] with an economy-size QR
   srf = size(rhs)[1]
   srs = size(rhs)[2]
   sr1s = size(rhs1)[2]
   rr = zeros(srf,srs+sr1s)
   rr[1:srf,1:srs] = rhs
   rr[1:srf,srs+1:srs+sr1s] = rhs1
   # Julia qr decomposition is always "economy size"
   U,beta=qr(rr)


   ibeta=inv(beta[1:s,1:s]);
   beta = beta[1:sh,1:sh];
   beta2=beta*beta';
   H=zeros((m+1)*s,m*s);
   T=zeros((m+1)*s,m*s);
   L=zeros((m+1)*s,m*s);
#   println("      it        backward err\n")
   global rho, js, j
   for j=1:m
     jms=(j-1)*s+1
     j1s=(j+1)*s
     js=j*s
     js1=js+1
     jsh=(j-1)*s+sh

     # Expand the basis
     # multiply by A
     Up = zeros(n,s)
     Up[1:n,1:sh] = LE\(A*(LE'\U[:,jms:jsh]))
     # solve with A

     Up[1:n,sh+1:s] = LE'*(UA\(LA\(LE*U[1:n,jsh+1:js])))

     # orthogonalize the new basis block wrt all the previous ones by modified gram
     for l=1:2
        k_min=max(1,j-k_max);
        for kk=k_min:j
            k1=(kk-1)*s+1
            k2=kk*s
            coef= U[1:n,k1:k2]'*Up
            H[k1:k2,jms:js] = H[k1:k2,jms:js]+ coef
            Up = Up - U[:,k1:k2]*coef
        end
      end

      if (j<=m)
        Up,H[js1:j1s,jms:js] = qr(Up);
        hinv=inv(H[js1:j1s,jms:js]);
      end


      ###############################################################
      # Recover the columns of T=U'*A*U (projection of A onto the space) from
      # the colums of H.
      # REMARK: we need T as coefficient matrix of the projected problem.
      I=eye(js+s)

      if (j==1)
        HL = zeros(s+sh,2*sh)
        HL[1:s+sh,1:sh] = H[1:s+sh,1:sh]/ibeta[1:sh,1:sh]
        HL[1:s+sh,sh+1:2*sh] = eye(s+sh,sh)/ibeta[1:sh,1:sh]
        L[1:j*s+sh,(j-1)*sh+1:j*sh] = HL*ibeta[1:s,sh+1:s];
      else
        L[1:j*s+s,(j-1)*sh+1:j*sh] = L[1:j*s+s,(j-1)*sh+1:j*sh] + H[1:j*s+s,jms:jms-1+sh]*rho;
      end

      odds = [odds; jms:(jms-1+sh)]   # store the odd block columns
      evens = 1:js
      flag = trues(size(evens))
      flag[odds] = false
      evens = evens[flag]
      T[1:js+s,odds]=H[1:js+s,odds]   #odd columns

      T[1:js+sh,evens]=L[1:js+sh,1:j*sh]   #even columns
      L[1:j*s+s,j*sh+1:(j+1)*sh] = ( I[1:j*s+s,(js-sh+1):js]- T[1:js+s,1:js]*H[1:js,js-sh+1:js])*hinv[sh+1:s,sh+1:s]
      rho = hinv[1:sh,1:sh]\hinv[1:sh,sh+1:s]

      #################################################################

      # Solve the projected problem by Bartels-Stewart
      # Do "type lyap" from command window if interested
      k=j

      Y = lyap((T[1:js,1:js]),eye(k*s,sh)*beta2*eye(k*s,sh)')

      # safeguard to preserve symmetry
      Y = (Y+Y')/2

      # Compute the residual norm. See the article by Valeria
      cc = zeros(j1s-js1+1,s)
      cc[:,1:s-sh]=H[js1:j1s,js-s+1:js-sh]
      cc[:,s-sh+1:s]=L[js1:j1s,(j-1)*sh+1:j*sh]

      nrmx = vecnorm(Y)

      er2=[er2;sqrt2*vecnorm(cc*Y[js-s+1:js,:])/(nrmb+singE*nrma*nrmx)]

#      @printf("It: #d, Current relative residual norm: #10.5e \n",k,er2[k])

      if (er2[k]<tol)
        break
      else
        su = size(U)[2]
        sup = size(Up)[2]
        newU = zeros(n,su+sup)
        newU[1:n,1:su]=U
        newU[1:n,su+1:su+sup]=Up
        U = newU
      end
    end
    # Done
    # reduce solution rank if needed
    sY,uY=eig(Y)
    sY=sort(sY)
    id=sortperm(sY)
    sY=flipdim(sY,1)
    uY=uY[:,id[end:-1:1]]
    is = 0
    for ii in 1:size(sY)[1]
      if abs(sY[ii])>tolY
        is = is+1
      end
    end

    Y0 = uY[:,1:is]*diagm(sqrt(sY[1:is]))
    Z = LE'\U[1:n,1:js]*Y0

    er2=er2[1:j]

    total_time=toq()

#    println("   its           comp.res.   space dim.   CPU Time\n")
#    @printf("#10.5e  #10.5e   #10.5e  #10.5e \n",j,er2[j],js,total_time)
    #println("************* \n AT CONVERGENCE \n")
    #@printf("Its: #d, Computed res: #10.5e, Space dim: #d, CPU Time: #10.5e\n",j,er2[j],js,total_time)
    return Z, er2

end

end
