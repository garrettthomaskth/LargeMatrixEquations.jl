module matrixEqs

export kpik, rksm

typealias ScalarOrArray{T} Union{T, Array{T}}

###############################################
#
# kpik function definition
#
###############################################

function kpik(A::ScalarOrArray,B::ScalarOrArray,E=1;LE=1,m::Number=100,tol::Number=1e-9,tolY::Number=1e-12)
  # Julia code for K-PIK (Krylov-plus-inverted-Krylov)
  # Based on kpik.m avalible from V. Simoncini's website
  #
  # Approximately solve
  #
  #       A X E + E X A' + BB' = 0
  #
  # by means of the extended Krylov subspace method
  # Input
  #  A   coeff matrix, A < 0
  #  B   factor of rhs,   nxk matrix with k << n
  #  NAMED ARGUMENTS
  #  E   coeff matrix, spd, Defult: 1
  #  LE  lower triang factor of coeff matrix, Defult: 1
  #  *Note: This is an optional argument, if not provided it will be set to
  #  cholfact(E)[:L]*
  #  m   max space dimension, Defult: 100
  #  tol stopping tolerance, with stopping criterion
  #          ||LE\A X LE  + LE' X A'/LE'-LE\BB'/LE'||
  #          ----------------------------------------  < tol
  #      ||LE\BB'/LE'|| + ||E^{-1}|| ||A|| ||LE'X LE ||
  #      computed in a cheap manner. Defult: 1e-9
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
  # * To solve
  #
  #       A X + X A' + BB' = 0
  #
  # Use kpik(A,B) as E is set to 1 by Defult
  #
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

  @assert(isdefined(:vecnorm),"Your julia version is too old. vecnorm() not defined")

  tic()

  #### Check if we solve the general lyapunov equation and user did not provide LE
  if E != 1 && LE == 1
    #### If this is the case, calculate LE
    LE = cholfact(E)[:L]
  end

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
    singE=1
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
   rr = [ rhs rhs1 ]

   # Julia qr decomposition is always "economy size"
   U,beta=qr(rr)
   U = U[1:n,1:s]

   ibeta=inv(beta[1:s,1:s]);
   beta = beta[1:sh,1:sh];
   beta2=beta*beta';
   H=zeros((m+1)*s,m*s);
   T=zeros((m+1)*s,m*s);
   L=zeros((m+1)*s,m*s);
   println("      it        backward err\n")
   local js, j, rho
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
      Iden=eye(js+s)

      if (j==1)
        L[1:j*s+sh,(j-1)*sh+1:j*sh] = [H[1:s+sh,1:sh]/ibeta[1:sh,1:sh] eye(s+sh,sh)/ibeta[1:sh,1:sh]]*ibeta[1:s,sh+1:s];
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
      L[1:j*s+s,j*sh+1:(j+1)*sh] = ( Iden[1:j*s+s,(js-sh+1):js]- T[1:js+s,1:js]*H[1:js,js-sh+1:js])*hinv[sh+1:s,sh+1:s]
      rho = hinv[1:sh,1:sh]\hinv[1:sh,sh+1:s]

      #################################################################

      # Solve the projected problem by Bartels-Stewart
      # Do "type lyap" from command window if interested

      Y = lyap((T[1:js,1:js]),eye(j*s,sh)*beta2*eye(j*s,sh)')

      # safeguard to preserve symmetry
      Y = (Y+Y')/2

      # Compute the residual norm. See the article by Valeria

      cc = [H[js1:j1s,js-s+1:js-sh] L[js1:j1s,(j-1)*sh+1:j*sh]]

      nrmx = vecnorm(Y)

      er2=[er2;sqrt2*vecnorm(cc*Y[js-s+1:js,:])/(nrmb+singE*nrma*nrmx)]

      @printf("It: %d, Current relative residual norm: %10.5e \n",j,er2[j])

      (er2[j]<tol) ? break : U = [U Up]

    end
    # Done
    # reduce solution rank if needed
    sY,uY=eig(Y)
    id=sortperm(sY)
    sY=sort(sY)

    sY=flipdim(sY,1)
    uY=uY[:,id[end:-1:1]]
    is = 0
    for ii in 1:size(sY)[1]
      if abs(sY[ii])>tolY
        is = is+1
      end
    end

    Y0 = uY[:,1:is]*diagm(sqrt(sY[1:is]))
    Z = LE'\(U[1:n,1:js]*Y0)
    total_time=toq()
    er2=er2[1:j]
    println("   its           comp.res.   space dim.   CPU Time\n")
    @printf("%10.5e  %10.5e   %10.5e  %10.5e \n",j,er2[j],js,total_time)
    return Z, er2

end



###############################################
#
# rksm function definition
#
###############################################


function rksm(A::ScalarOrArray,B::ScalarOrArray,E::ScalarOrArray=1;EL::ScalarOrArray=1,
  s1::Number=NaN,emax::Number=NaN,m::Int=100,tol::Number=1e-9,ch::Bool=true,tolY::Number=1e-12)
# Based on rksm.m on Valeria Simoncini's website
#
# Approximately Solve
#                A X E + E X A' + BB' = 0
#
# by the Rational Krylov subspace method
# (Galerkin condition onto the Rational Krylov subspace)
# This code performs system solves with (A-s E)
#
# Input:
#
# A, E  coeff. matrices. A<0,  E is spd
# EL   Cholesky lower factor of E
# B     rhs factor
# m       max space dimension allowed
# tol     stopping tolerance, with stopping criterion
#          ||LE\A X LE  + LE' X A'/LE'-LE\BB'/LE'||
#          ----------------------------------------  < tol
#      ||LE\BB'/LE'|| + ||E^{-1}|| ||A|| ||LE'X LE ||
#         computed in a cheap manner
# s1,smax estimates for real spectral interval
#         associated with field of values of (A,E)
# ch      ch=true  complex poles  ch=false real poles
# tolY    truncation tolerance for final solution, e.g., tolY=1e-12
#
# Output:
#
# Z    factor of approximate solution  X = Z Z'
#      If ch=true,  Z may be complex
# nrmrestot  history of residuals
#
# Hints:
# 1) Before the call permute entries of A, E and B so as to
#    limit fill-in in the system solves
# 2) Provide "comfortable" (loose bounds) estimates s1, emax
#
#  When using this code, please cite the following reference:
#
#  V. Druskin and V. Simoncini,
#  Adaptive rational Krylov subspaces for large-scale dynamical systems
#  Systems & Control Letters, 60 (2011), pp. 546-560.
#
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
#FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
#COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
#IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
@assert(isdefined(:vecnorm),"Your julia version is too old. vecnorm() not defined")

tic()

symm = norm(A-E-(A-E)',1) < 1e-14
# If not symmetric, we may have complex poles.
# Therefore we must use complex types
symm ? typ = Float64 : typ = Complex{Float64}

# If user does not give smallest generalized eigenvalue, calculate it
if (isnan(s1))
  (symm) && (s1 = eigmin(inv(E)*(-A)))
  (!symm) && (s1 = minimum(abs(eigvals(inv(E)*(-A)))))
end

# If user does not give largest generalized eigenvalue, calculate it
if (isnan(emax))
  (symm) && (emax = eigmax(inv(E)*(-A)))
  (!symm) && (emax = maximum(abs(eigvals(inv(E)*(-A)))))
end

if (E != 1 && EL == 1)
  EL = cholfact(E)[:L]
end

n=size(A,1)
B=full(B)
p=size(B,2)
Iden=speye(p)
O=0*Iden
uno=ones(1,p)
Lres = convert(Array{typ,2},EL\B)

V,irr=qr(Lres)
rr=inv(irr)
nrmb=vecnorm(inv(rr))^2
beta=V'*Lres
beta2=beta*beta'
s=typ[]
print("     no its     backward error\n")
VV=Array(typ,n,p*(m+2))
VV[1:n,1:p]=V
H=Array(typ,p*(m+2),p*(m+1))
nrmrestot=[]
nrma=vecnorm(A)


if (vecnorm(E-speye(n))>1e-14)
  condestE=cond(E);
  singE=condestE/vecnorm(E);
else
  singE=1;
end


newAv=EL\(A*(EL'\V));

K=V'*newAv;
push!(s,s1);
eH=eig(K)
eHpoints = sort([s1,emax])
snew=newpolei(eHpoints,eH[1],s1*uno',typ);

if real(snew)<0
   snew=-real(snew)+im*imag(snew);
end
push!(s,snew);

# additional steps
cmplxflag=false;

i=0;
local Y, j1s
while i < m
  local i1, j1, js
  i=i+1;

  paired=0;
  while (paired==0)

    i1=i+1;
    w=EL*V;
    E==1 ? wrk = (A-snew*diagm(ones(n),0))\w : wrk = (A-snew*E)\w;
    wrk= EL'*wrk;

    # Gram-Schmidt step
    jms=(i-1)*p+1;
    j1s=(i+1)*p;
    js=i*p;
    js1=js+1;
    for it=1:2
      for kk=1:i
        k1=(kk-1)*p+1;
        k2=kk*p;
        gamma=VV[1:n,k1:k2]'*wrk;
        H[k1:k2,jms:js] = H[k1:k2,jms:js]+ gamma;
        wrk = wrk - VV[:,k1:k2]*gamma;
      end
    end

    V, H[js1:j1s,jms:js]=qr(wrk);

    if (cmplxflag)
      snew=conj(snew);
      push!(s,snew)
      cmplxflag=false;
      newAv=EL\(A*(EL'\V));
      g = VV[1:n,1:js]'*newAv;
      g1 = g;
      g2 = V'*(EL\(A*(EL'\VV[1:n,1:js])));
      g3 = V'*(EL\(A*(EL'\V)));
      K = [K g1; g2 g3];
      VV[1:n,js+1:j1s]=V;
      i=i+1;
    else
      paired=1
    end
  end

  ih1=i1;
  ih=i;
  newAv=EL\(A*(EL'\V));
  g = VV[1:n,1:js]'*newAv;

  if (symm)
    K=(K+K')/2;
  end

  rhs2=speye(ih*p,p)*beta2*speye(ih*p,p)';
  Y = lyap(K,rhs2);
  nrmx = vecnorm(Y);

  # computed residual   (exact, in exact arithmetic)
  u1=newAv-VV[1:n,1:js]*g;
  d=-VV[1:n,1:js]*(Y*(H[1:ih*p,1:ih*p]'\[zeros(p*(ih-1),p);Iden])*H[p*ih+1:p*ih1,p*ih-p+1:p*ih]');
  U=[-V*s[end]  d u1 ];
  extra,rr=qr(full(U));
  nrmres=vecnorm(rr*sparse([O Iden O; Iden O Iden; O Iden O ])*rr')/(nrmb+singE*nrma*nrmx);
  nrmrestot=[nrmrestot; nrmres];

  println([i,nrmres])

  (nrmres<tol) && break

  # New poles and zeros
  eH=eig(K)[1];
  symm ? sort!(eH) : eH[sortperm(abs(eH))]

  eHorig=eH;

   if (ch)
     # Complex poles. Compute set for next complex pole of r_m

      if (countnz(imag(eH))>0 && length(eH)>2) # Roots lambdas come from convex hull too
        eH=full([eH;-emax]);
        eH=convhull(eH)
        ieH=length(eH);
        missing=ih*p-ieH;
        while missing>0                        # include enough points from the border
          neweH=(eH[1:ieH-1]+eH[2:ieH])/2;
          missing=ih*p-length(eH);
          eH=[eH;neweH];
        end
        eHpoints=-eH;
        eH=eHorig;
      else     # if all real eigs, no convex hull possible
        eHpoints = sort([s1; emax.';-real(eH)]);
      end


   else   # Real poles s from real set. Compute complex roots of r_m via Ritz convex hull
       if (countnz(imag(eH))>0 && length(eH)>2)    # Roots lambdas come from convex hull too
         eH=full([eH;-s1;-emax.']);
         eH=convhull(eH)
         ieH=length(eH);
         missing=ih*p-ieH;
         while missing>0 # include enough points from the border
           neweH=(eH[1:ieH-1]+eH[2:ieH])/2;
           eH=[eH;neweH];
           missing=ih*p-length(eH);
         end
         eH=eH[1:ih*p];
       end
        eHpoints = sort([s1; emax.';-real(eH)]);
        eH=eHorig;
   end

   gs=kron(s[2:i+1].',uno)';

   snew = newpolei(eHpoints,eH,gs,typ);
   (real(snew)<0) && (snew=-real(snew)+im*imag(snew)); #safeguard strategy

   # If pole is complex, include its conjugate
   (imag(snew) !=0) && (cmplxflag=true);

   push!(s,snew)

   g1 = g;
   g2 = V'*(EL\(A*(EL'\VV[1:n,1:js])));
   g3 = V'*(EL\(A*(EL'\V)));
   K = [K g1; g2 g3];
   VV[1:n,js+1:j1s]=V;

end;

# Done
# Reduce rank of solution, if needed
sY,uY=eig(Y)
symm ? id = sortperm(sY) : id = sortperm(abs(sY))

sY=sY[id]
sY=flipdim(sY,1)
uY=uY[:,id[end:-1:1]]
is = 0
for ii in 1:size(sY)[1]
  if abs(sY[ii])>tolY
    is = is+1
  end
end

Y0 = uY[:,1:is]*diagm(sqrt(sY[1:is]))
Z = EL'\(VV[:,1:size(Y0,1)]*Y0)

RKStotal_time=toq();


@printf("Space dim %d  Solution rank %d  time %10.5e \n",j1s,is,RKStotal_time);
return Z,nrmrestot

end



##################################
# Auxiliary routines

##################################
function ratfun(x,eH,s,typ)
r = Array(typ,length(x))

for j=1:length(x)
  xj = x[j]*ones(length(s))
  r[j]=abs(prod( (xj-s)./(xj-eH) ));
end
return r
end

##################################
function newpolei(eHpoints,eH,s,typ)
snew=Array(typ,length(eHpoints)-1)

for j=1:length(eHpoints)-1

    sval = linspace(real(eHpoints[j]),real(eHpoints[j+1]),200)
    (typ == Complex{Float64}) && (sval = sval + im*linspace(imag(eHpoints[j]),imag(eHpoints[j+1]),200))

    jx = indmax(abs(ratfun(sval,eH,s,typ)))

    snew[j]=sval[jx];
end
jx=indmax(abs(ratfun(snew,eH,s,typ)))
snew=snew[jx];
return snew
end



#######################################
function convhull(pnts)
    #Based on code from (https://github.com/intdxdt/convexhull.jl)

    # Function to compute 2-D convex hull
    T = eltype(pnts) #get point type
    N = length(pnts) #number of pnts
    # sort the points lexicographically.
    # copy points into mutable container

    #trivial case 0 or 1 point
    length(pnts) <= 1 && (return pnts)

    #sort points lexicographically
    sort!(pnts, lt=lt2d)

    #function to orient boundry using 2D cross product of OA and OB vectors,
    #i.e. z-component of their 3D cross product.
    #Returns a positive value, if OAB makes a counter-clockwise turn,
    #negative for clockwise turn, and zero if the points are collinear.
    orient(b, pnt) = (real(b[end]-b[end-1]) * imag(pnt - b[end-1])) -
                      (imag(b[end-1] - b[end-1]) * real(pnt - b[end-1]))

    #build lower hull
    lower = T[]
    for i in 1:N
        pnt = pnts[i]
        while length(lower) >= 2 && orient(lower, pnt) <= 0
            pop!(lower)
        end
        push!(lower, pnt)
    end

    #build upper hull
    upper = T[]
    for i in N:-1:1
        pnt = pnts[i]
        while length(upper) >= 2 && orient(upper, pnt) <= 0
            pop!(upper)
        end
        push!(upper, pnt)
    end

    # omit end points :
    #   - reapeated at the beginning of lower & upper
    length(upper) > 0 && pop!(upper)
    length(lower) > 0 && pop!(lower)

    # two reapeated points
    (length(upper)==1 && length(lower)==1) &&
    isequal(upper, lower) && pop!(upper)

    # concat lower and upper hull.
    hull = [lower; upper]

    #close ring
    length(hull) > 1 && push!(hull, hull[1])

    return hull
end

#less than comparator
function lt2d(a, b)
    dx, dy =  real(a - b), imag(a - b);
    dx != 0 &&  (return <(dx, 0))
    return <(dy, 0)
end



end
