# Translation into Julia of
# V. Simoncini
# A new iterative method for solving large-scale Lyapunov matrix equations,
# SIAM J.  Scient. Computing, v.29, n.3 (2007), pp. 1268-1288.
# Davide's Edit: doesn't manage the solution of a generalized Lyapunov equation
function kpikSimp(A,B,m,tol,tolY)
  tic()
  rhs=B
  nrmb=vecnorm(rhs)^2
  sqrt2=sqrt(2)
  er2=zeros(m,1)

  n,sh=size(rhs)
  Y=[]
  odds=[]
  er2=[]
  global rho, js, j
  # factorize A just once
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
   s=2*sh

   #Solve with A: rhs1=A^{-1}B
   rhs1=UA\(LA\rhs)

   # Orthogonalize [B,A^{-1}B] with an economy-size QR
   srf = size(rhs)[1]
   srs = size(rhs)[2]
   sr1s = size(rhs1)[2]
   rr = zeros(srf,srs+sr1s)
   rr[1:srf,1:srs] = rhs
   rr[1:srf,srs+1:srs+sr1s] = rhs1
   U,beta,p=qr(rr, Val{true})
   #U = U[1:n,1:s]




   ibeta=inv(beta[1:s,1:s])
   beta = beta[1:sh,1:sh]
   beta2=beta*beta'

   # Preallocate
   H=zeros((m+1)*s,m*s)
   T=zeros((m+1)*s,m*s)
   L=zeros((m+1)*s,m*s)

   for j=1:m
     jms=(j-1)*s+1
     j1s=(j+1)*s
     js=j*s
     js1=js+1
     jsh=(j-1)*s+sh

     # Expand the basis
     # multiply by A
     Up = zeros(n,s)
     Up[1:n,1:sh] = A*U[:,jms:jsh]
     # solve with A

     Up[1:n,sh+1:s] = UA\(LA\U[1:n,jsh+1:js])


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

      # now the new basis block is orthogonal wrt the previous ones, but its
      # columns are not orthogonal wrt each other --> economy-size QR
      #print(Up)
      #print("\n \n")
      Up,H[js1:j1s,jms:js] = qr(Up)
      #print(Up)
      #print("\n \n")
      hinv=inv(H[js1:j1s,jms:js])

      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # Recover the columns of T=U'*A*U (projection of A onto the space) from
      # the colums of H.
      # REMARK: we need T as coefficient matrix of the projected problem.
      I=eye(js+s)

      if (j==1)
        HL = zeros(s+sh,2*sh)
        HL[1:s+sh,1:sh] = H[1:s+sh,1:sh]/ibeta[1:sh,1:sh]
        HL[1:s+sh,sh+1:2*sh] = speye(s+sh,sh)/ibeta[1:sh,1:sh]
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

      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      # Solve the projected problem by Bartels-Stewart
      # Do "type lyap" from command window if interested
      k=j
      #print("\n")
      #print(T[1:js,1:js])
      #print("\n")
      #print(eye(k*s,sh)*beta2*eye(k*s,sh)')
      #print("\n")
      Y = lyap((T[1:js,1:js]),eye(k*s,sh)*beta2*eye(k*s,sh)')

      # safeguard to preserve symmetry
      Y = (Y+Y')/2

      # Compute the residual norm. See the article by Valeria
      cc = zeros(j1s-js1+1,s)
      cc[:,1:s-sh]=H[js1:j1s,js-s+1:js-sh]
      cc[:,s-sh+1:s]=L[js1:j1s,(j-1)*sh+1:j*sh]
      #HL = [H[js1:j1s,js-s+1:js-sh]; L[js1:j1s,(j-1)*sh+1:j*sh]]
      #d1 = j1s-js1+1
      #d2 = convert(Int64,size(HL)[1]/d1)

      er2=[er2;sqrt2*vecnorm(cc*Y[js-s+1:js,:])/nrmb]

      @printf("It: %d, Current relative residual norm: %10.5e \n",k,er2[k])

      if (er2[k]<tol)
        break
      else
        su = size(U)[2]
        sup = size(Up)[2]
        newU = zeros(n,su+sup)
        newU[1:n,1:su]=U
        newU[1:n,su+1:su+sup]=Up
        U = newU
        #print(U)
        #print("\n \n")
      end
    end
    # Done
    # reduce solution rank if needed
    sY,uY=eig(Y)
    # !!!!!!!!!!! sort returns different than matlab
    sY=sort(sY)
    id=sortperm(sY)
    sY=flipud(sY)
    uY=uY[:,id[end:-1:1]]
    is = 0
    for ii in 1:size(sY)[1]
      if abs(sY[ii])>tolY
        is = is+1
      end
    end

    Y0 = uY[:,1:is]*diagm(sqrt(sY[1:is]))
    Z = U[1:n,1:js]*Y0

    er2=er2[1:j]

    total_time=toq()
    println("************* \n AT CONVERGENCE \n")
    @printf("Its: %d, Computed res: %10.5e, Space dim: %d, CPU Time: %10.5e\n",j,er2[j],js,total_time)
end
