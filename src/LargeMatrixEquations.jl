module LargeMatrixEquations

export kpik, rksm, lp_lradi, lp_para


"""
    `kpik(A,B,E=1;<keyword arguments>)`

Julia code for K-PIK (Krylov-plus-inverted-Krylov)

Based on kpik.m avalible from V. Simoncini's website
(http://www.dm.unibo.it/~simoncin/software.html)
and includes many of the same comments and much of the same description.

Approximately solve

     A X E + E X A' + BB' = 0

by means of the extended Krylov subspace method

ARGUMENTS

    'A' : coeff matrix, A < 0

    'B' : factor of rhs,   nxk matrix with k << n

NAMED ARGUMENTS

    'E' : coeff matrix, spd, Defult: 1

    'LE' : lower triang factor of coeff matrix, Defult: 1 Note: This is an
    optional argument, if not provided it will be set to cholfact(E)[:L]

    'm' : max space dimension, Defult: 100

    'tol' : stopping tolerance based on the backwards error,
    with stopping criterion

      ||LE\\A X LE  + LE' X A'/LE'-LE\\BB'/LE'||
    -----------------------------------------------   < tol
    ||LE\\BB'/LE'|| + ||E^{-1}|| ||A|| ||LE'X LE ||

    computed in a cheap manner. Defult: 1e-9

Output

    'Z' : solution factor s.t.  X = Z Z'

    'er2' : history of scaled residual, as above

Comments

    1. The projected solution is computed at each iteration
    As an alternative, a periodic computation could be considered.

    2. This code performs a factorization of A. As an alternative,
    iterative solves could be considered.

    3. To solve

     A X + X A' + BB' = 0

     Use kpik(A,B) as E is set to 1 by Defult

If you use this code, please cite the following article:

    V. Simoncini
    A new iterative method for solving large-scale Lyapunov matrix equations,
    SIAM J.  Scient. Computing, v.29, n.3 (2007), pp. 1268-1288.
"""
function kpik(A,B,E=1;LE=1,m=100,tol=1e-9,tolY=1e-12,infoV=true)

  @assert(isdefined(:vecnorm),"Your julia version is too old. vecnorm() not defined")

  infoV && tic()

  #### Check if we solve the general lyapunov equation and user did not provide LE
  if E != 1 && LE == 1
    #### If this is the case, calculate LE
    if issparse(E)
      LE = try
            sparse(cholfact(E,perm=1:size(E,1))[:L])
           catch
            error("E must be SPD")
          end
    else
      LE = try
            cholfact(E)[:L]
          catch
            error("E must be SPD")
          end
    end
  end

  const rhs=LE\B
  const nrmb=vecnorm(rhs)^2
  const nrma=vecnorm(A)
  const sqrt2=sqrt(2)
  er2=zeros(m,1)

  const n,sh=size(rhs)

  Y=[]
  odds=[]
  er2=[]

  if E == 1
    const condestE=1
    const singE=condestE/vecnorm(eye(n))
  elseif (vecnorm(E-speye(n))>1e-14)
    const condestE=cond(full(E))
    const singE=condestE/vecnorm(E)
  else
    const singE=1
  end

  if norm(A-A',1)<1e-14
     if issparse(A)
       cfA = cholfact(-A,perm=1:size(A,1))
       const LA = sparse(cfA[:L])
       const UA = -LA'
     else
       const UA = chol(-A)
       const LA = -UA'
     end
     infoV && println("A sym. Completed Chol factorization\n")
     const k_max =2
   else
     luA=lufact(full(A),Val{false})
     const LA = luA[:L]
     const UA = luA[:U]


     infoV && println("A nonsym. Completed LU factorization\n")
     const k_max = m
   end

   const s=2*sh
   #rhs1=LE'*(UA\(LA\(LE*rhs)))

   # Orthogonalize [B,A^{-1}B] with an economy-size QR
   rr = [ rhs LE'*(UA\(LA\(LE*rhs))) ]

   # Julia qr decomposition is always "economy size"
   U,beta=qr(rr)
   U = U[1:n,1:s]

   const ibeta=inv(beta[1:s,1:s])
   beta = beta[1:sh,1:sh]
   const beta2=beta*beta'
   H=zeros((m+1)*s,m*s)
   T=zeros((m+1)*s,m*s)
   L=zeros((m+1)*s,m*s)
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
        k_min=max(1,j-k_max)
        for kk=k_min:j
            k1=(kk-1)*s+1
            k2=kk*s
            coef= U[1:n,k1:k2]'*Up
            H[k1:k2,jms:js] = H[k1:k2,jms:js]+ coef
            Up = Up - U[:,k1:k2]*coef
        end
      end

      if j<=m
        Up,H[js1:j1s,jms:js] = qr(Up)
        hinv=inv(H[js1:j1s,jms:js])
      end


      ###############################################################
      # Recover the columns of T=U'*A*U (projection of A onto the space) from
      # the colums of H.
      # REMARK: we need T as coefficient matrix of the projected problem.
      Iden=speye(js+s)

      if (j==1)
        L[1:j*s+sh,(j-1)*sh+1:j*sh] = [H[1:s+sh,1:sh]/ibeta[1:sh,1:sh] eye(s+sh,sh)/ibeta[1:sh,1:sh]]*ibeta[1:s,sh+1:s]
      else
        L[1:j*s+s,(j-1)*sh+1:j*sh] = L[1:j*s+s,(j-1)*sh+1:j*sh] + H[1:j*s+s,jms:jms-1+sh]*rho
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

      infoV && println("KPIK It: $j -- Current Backwards Error: $(er2[j])")

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
    er2=er2[1:j]

    if infoV
      println("its  Back. Error            space dim. CPU Time")
      println("$j    $(er2[j])  $js          $(toq())")
    end

    return Z, er2

end


"""
  `rksm(A,B,E=1;LE=1,s1=NaN,emax=NaN,m=100,tol=1e-9,ch=true,tolY=1e-12,infoV=true)`

Based on rksm.m on Valeria Simoncini's website
(http://www.dm.unibo.it/~simoncin/software.html)
and includes many of the same comments and much of the same description

Approximately Solve

    A X E + E X A' + BB' = 0

by the Rational Krylov subspace method
(Galerkin condition onto the Rational Krylov subspace)
This code performs system solves with (A-s E)

Input:

    A, E  coeff. matrices. A<0,  E is spd

    LE   Cholesky lower factor of E

    B     rhs factor

    m       max space dimension allowed

    tol     stopping tolerance based on the backwards error,
            with stopping criterion

                ||LE\A X LE  + LE' X A'/LE'-LE\BB'/LE'||
            ----------------------------------------------  < tol
            ||LE\BB'/LE'|| + ||E^{-1}|| ||A|| ||LE'X LE ||

            computed in a cheap manner

    s1,smax   estimates for real spectral interval associated with
              field of values of (A,E)

    ch      ch=true  complex poles  ch=false real poles

    tolY    truncation tolerance for final solution, e.g., tolY=1e-12

    infoV   If true, print out backwards error every iteration and a summary
            once the program finishes.

Output:

    Z    factor of approximate solution  X = Z Z'
         If ch=true,  Z may be complex

    nrmrestot  history of residuals

Hints:

    1. Before the call permute entries of A, E and B so as to
    limit fill-in in the system solves

    2. Provide "comfortable" (loose bounds) estimates s1, emax

When using this code, please cite the following reference:

    V. Druskin and V. Simoncini,
    Adaptive rational Krylov subspaces for large-scale dynamical systems
    Systems & Control Letters, 60 (2011), pp. 546-560.


#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
#FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
#COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
#IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
function rksm(A,B,E=1;LE=1,s1=NaN,emax=NaN,m=100,tol=1e-9,ch=false,tolY=1e-12,
                                                          infoV=true)

@assert(isdefined(:vecnorm),"Your julia version is too old. vecnorm() not defined")

infoV && tic()

const symm = norm(A-E-(A-E)',1) < 1e-14
# If not symmetric, we may have complex poles.
# Therefore we must use complex types
symm ? (const typ = Float64) : (const typ = Complex{Float64})

const n=size(A,1)
const p=size(B,2)

E==1 ? ee = eye(n) : ee = E
# If user does not give smallest generalized eigenvalue, calculate it
if (isnan(s1))
  (symm) && (s1 = eigmax(inv(full(ee))*full(-A)))
  (!symm) && (s1 = maximum(abs(eigvals(inv(full(ee))*full(-A)))))
end

# If user does not give largest generalized eigenvalue, calculate it
if (isnan(emax))
  (symm) && (emax = eigmin(inv(full(ee))*full(-A)))
  (!symm) && (emax = minimum(abs(eigvals(inv(full(ee))*full(-A)))))
end

if (E != 1 && LE == 1)
  if issparse(E)
    LE = try
          LE = sparse(cholfact(E,perm=1:size(E,1))[:L])
        catch
          error("E must be Symmetric Positive Definite")
        end
  else
    LE = try
          cholfact(E)[:L]
        catch
          error("E must be Symmetric Positive Definite")
        end
 end
end

const Iden=eye(p)
const O=0*Iden
const uno=ones(1,p)
const Lres = convert(Array{typ,2},LE\B)

V,irr=qr(Lres)
rr=inv(irr)
const nrmb=vecnorm(inv(rr))^2
const beta=V'*Lres
const beta2=beta*beta'
s=typ[]

VV=zeros(n,p*(m+2))
(typ == Float64) && (VV=VV*1.)
(typ == Complex{Float64}) && (VV=VV*1.*im)
VV[1:n,1:p]=V

H=zeros(p*(m+2),p*(m+1))
(typ == Float64) && (H=H*1.)
(typ == Complex{Float64}) && (H=H*1.*im)

nrmrestot=[]
const nrma=vecnorm(A)


if E == 1
  const condestE=1
  const singE=condestE/vecnorm(eye(n))
elseif (vecnorm(E-speye(n))>1e-14)
  const condestE=cond(full(E))
  const singE=condestE/vecnorm(full(E))
else
  const singE=1
end

newAv=LE\(A*(LE'\V))

K=V'*newAv
push!(s,s1)
eH=eig(K)
eHpoints = sort([s1,emax])
snew=newpolei(eHpoints,eH[1],s1*uno',typ)

if real(snew)<0
   snew=-real(snew)+im*imag(snew)
end
push!(s,snew)

# additional steps
cmplxflag=false

i=0


local Y, j1s
while i < m
  local i1, j1, js
  i=i+1

  paired=0
  while (paired==0)

    i1=i+1
    w=LE*V
    #E==1 ? wrk = (A-snew*diagm(ones(n),0))\w : wrk = (A-snew*E)\w
    wrk = (A-snew*speye(n)*E)\w
    wrk= LE'*wrk

    # Gram-Schmidt step
    jms=(i-1)*p+1
    j1s=(i+1)*p
    js=i*p
    js1=js+1
    for it=1:2
      for kk=1:i
        k1=(kk-1)*p+1
        k2=kk*p
        gamma=VV[1:n,k1:k2]'*wrk
        H[k1:k2,jms:js] = H[k1:k2,jms:js]+ gamma
        wrk = wrk - VV[:,k1:k2]*gamma
      end
    end
    V, H[js1:j1s,jms:js]=qr(wrk)
    if (cmplxflag)
      snew=conj(snew)
      push!(s,snew)
      cmplxflag=false
      newAv=LE\(A*(LE'\V))
      g = VV[1:n,1:js]'*newAv
      g1 = g
      g2 = V'*(LE\(A*(LE'\VV[1:n,1:js])))
      g3 = V'*(LE\(A*(LE'\V)))
      K = [K g1; g2 g3]
      VV[1:n,js+1:j1s]=V
      i=i+1
    else
      paired=1
    end
  end

  ih1=i1
  ih=i
  newAv=LE\(A*(LE'\V))
  g = VV[1:n,1:js]'*newAv

  if (symm)
    K=(K+K')/2
  end

  rhs2=speye(ih*p,p)*beta2*speye(ih*p,p)'
  Y = lyap(K,rhs2)
  nrmx = vecnorm(Y)

  # computed residual   (exact, in exact arithmetic)
  u1=newAv-VV[1:n,1:js]*g
  d=-VV[1:n,1:js]*(Y*(H[1:ih*p,1:ih*p]'\[zeros(p*(ih-1),p);Iden])*H[p*ih+1:p*ih1,p*ih-p+1:p*ih]')

  U=[-V*s[end]  d u1 ]
  extra,rr=qr(full(U))

  nrmres=vecnorm(rr*([O Iden O; Iden O Iden; O Iden O ])*rr')/(nrmb+singE*nrma*nrmx)
  push!(nrmrestot, nrmres)

  infoV && println("RKSM It: $i -- Current backwards error: $nrmres")
  (nrmres<tol) && break
  # New poles and zeros
  eH=eig(K)[1]
  symm ? sort!(eH) : eH[sortperm(abs(eH))]

  eHorig=eH

   if (ch)
     # Complex poles. Compute set for next complex pole of r_m

      if (countnz(imag(eH))>0 && length(eH)>2) # Roots lambdas come from convex hull too
        eH=full([eH;-emax])
        eH=convhull(eH)
        ieH=length(eH)
        missing=ih*p-ieH
        while missing>0                        # include enough points from the border
          neweH=(eH[1:ieH-1]+eH[2:ieH])/2
          missing=ih*p-length(eH)
          eH=[eH;neweH]
        end
        eHpoints=-eH
        eH=eHorig
      else     # if all real eigs, no convex hull possible
        eHpoints = sort([s1; emax.';-real(eH)])
      end


   else   # Real poles s from real set. Compute complex roots of r_m via Ritz convex hull
       if (countnz(imag(eH))>0 && length(eH)>2)    # Roots lambdas come from convex hull too
         eH=full([eH;-s1;-emax.'])
         eH=convhull(eH)
         ieH=length(eH)
         missing=ih*p-ieH
         while missing>0 # include enough points from the border
           neweH=(eH[1:ieH-1]+eH[2:ieH])/2
           eH=[eH;neweH]
           missing=ih*p-length(eH)
         end
         eH=eH[1:ih*p]
       end
        eHpoints = sort([s1; emax.';-real(eH)])
        eH=eHorig
   end

   gs=kron(s[2:i+1].',uno)'

   snew = newpolei(eHpoints,eH,gs,typ)
   (real(snew)<0) && (snew=-real(snew)+im*imag(snew)); #safeguard strategy

   # If pole is complex, include its conjugate
   (imag(snew) !=0) && (cmplxflag=true)

   push!(s,snew)

   g1 = g
   g2 = V'*(LE\(A*(LE'\VV[1:n,1:js])))
   g3 = V'*(LE\(A*(LE'\V)))
   K = [K g1; g2 g3]
   VV[1:n,js+1:j1s]=V

end


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
Z = LE'\(VV[:,1:size(Y0,1)]*Y0)

if infoV
  RKStotal_time=toq()
  println("Space dim $j1s  Solution rank $is  time $RKStotal_time")
end

return Z,nrmrestot

end



#################################################
#
#  ratfun and newpolei function definitions. These
#  are auxiliary routines for rksm based on Simoncini's
#  corresponding routines.
#
##################################################


function ratfun(x,eH,s,typ)
r = Array(typ,length(x))

for j=1:length(x)
  xj = x[j]*ones(length(s))
  r[j]=abs(prod( (xj-s)./(xj-eH) ))
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

    snew[j]=sval[jx]
end
jx=indmax(abs(ratfun(snew,eH,s,typ)))
snew=snew[jx]
return snew
end


"""
  `convexhull(pnts)`
Based on code from (https://github.com/intdxdt/convexhull.jl)
Function to compute 2-D convex hull
"""
function convhull(pnts)
    const T = eltype(pnts) #get point type
    const N = length(pnts) #number of pnts
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
    if (length(upper)==1 && length(lower)==1 && isequal(upper, lower))
      pop!(upper)
    end

    # concat lower and upper hull.
    hull =[lower; upper]

    #close ring
    length(hull) > 1 && push!(hull, hull[1])

    return hull
end

#less than comparator
function lt2d(a, b)
    dx, dy =  real(a - b), imag(a - b)
    dx != 0 &&  (return <(dx, 0))
    return <(dy, 0)
end




############################################################
#
#     Please note: The following functions lp_lradi, lp_nrmu, lp_para,
#                  lp_arn, lp_mnmx are based on the functions
#                  with the same name in LYAPACK 1.0 (lp_arn
#                  is based on lp_arn_m.m and lp_arn_p.m). The
#                  structure of these files is closely followed
#                  and even includes many of the same comments.
#
#############################################################


"""
Based on lp_lradi.m from LYAPACK 1.0

Low Rank Cholesky Factor Alternating Direction Implicit (LRCF-ADI) Method for
solving the stable Lyapunov equation:

  For tp = :B

    F * X + X * F' = -B * B',

  for tp = :C

    F' * X + X * F = -B' * B,

where F = A-Bf*Kf'.

The routine works in two modes (depending on the choice of zk)

    For zk = :Z, this routine deliveres a low rank factor Z,
    such that Z * Z' approximates X

    for zk = :K, it only generates the product K_out = Z * Z'* K
    without forming Z itself.

Input:

    A         matrix A
    B         matrix B (n-x-m if tp=:B or m-x-n if tp=:C,
              where  m << n !)
    p         l-vector with ADI shift parameters. Complex parameters must
              appear as conjugate complex pairs in consecutive order.
              If more than length(p) iterations need to be done, the
              parameters p[i] are reused in a cyclic manner. The user can calculate
              these values using lp_para.
    Bf        "feedback" matrix Bf
              Named argument with defult value [], so if not provided by user,
              the code solves
              A * X + X * A' = -B * B' or A' * X + X * A = -B'* B
              depending on choice of tp.
    Kf        "feedback" matrix Kf
              amed argument with defult value [], so if not provided by user,
              the code solves
              A * X + X * A' = -B * B' or A' * X + X * A = -B' * B
              depending on choice of tp.
    tp        (= :B or :C) named argument that determines the type of Lyapunov equation. Defult value is :B
    rc        (= :R or :C) named argument that determines whether the low rank factor Z
              must be real (:R) or if complex factors are allowed (:C)
              If p contains complex parameters, then the low rank factor
              Z is complex, too, although Z * Z' is real. A real low rank
              factor is determined from the complex data, if rc = :R.
              However, this requires some additional computation.
              If zk = :K, rc is ignored. Defult value is :C.
    K         named argument n-x-r matrix (where r should be small: r << n !).
              Defult value is [].
    max_it    stopping criterion: maximal number of LRCP-ADI steps
              Set to [] or +Inf (default) to avoid this criterion. Named argument
              With defult value of 100.
    min_res   stopping criterion: minimal relative residual norm. The
              iteration is stopped when res(i+1) <= min_res (See Remarks).
              Set min_res = [] or 0 (default) to avoid this criterion. Note:
              If min_res<=0 and with_rs=:N, the (often expensive) calculation
              of the residual norm is avoided, but, of course, res is not
              provided on exit. Named argument with defult value of 1e-6.
    with_rs   (= :S or :N) if with_rs = "S", the iteration is stopped,
              when the routine detects a stagnation of the residual norms,
              which is most likely the case, when roundoff errors rather
              than the approximation error start to dominante the residual
              norm. This implies that the residual norms are computed (which
              can be expensive). This criterion works quite well in practice,
              but is not absolutely sure. Use with_rs = "S" only if you
              want to compute the Lyapunov solution as accurate as possible
              (for a given machine precision).
              If with_rs = :N, this criterion is not used. Named argument with
              defult value :S.
    min_in    stopping criterion: This value limits the "minimal increase"
              in the matrix Z by the "new" columns . The iteration is
              terminated as soon as

                || Z_nc ||F^2 < min_in * || Z_new ||F^2

              holds for a certain number of consecutive iteration steps.
              Here, Z_nc are the currently computed "new" colums, which
              appended to the old iterate Z_old deliver the new iterate
              Z_new  = [ Z_old  Z_nc ]. Set min_in = 0 to avoid it.
              Default value is eps, the machine precision. min_in = []
              has the same effect. Note that this is an "adaptive"
              stopping criterion which does not require the
              (often expensive) computation of the residual norm. Named argument
              with defult value 0.
    infoV      Bool: Return infoVrmation given during the
              iteration. Named argument with default value true.

Output:

    Z         Z * Z' is a low rank approximation of X
              (Note that Z can be complex if rc=:C!)
    flag      the criterion which had stopped the iteration:
               = "I": max_it,
               = "R": min_res,
               = "S": with_rs,
               = "N": min_in,
    res       the relative residual norms attained in the course of
              the iterations (res[i+1] is the norm after the i-th step
              of the iteration!). See note in min_res.

Remarks:

    1. The eigenvalues of F must have negative real parts.

    2. The values in res correspond to the following "relative" norms

      tp = :B

        res(i+1) = ||F * Z_i * Z_i' + Z_i * Z_i' * F' + B * B'||F/|| B * B'||F

      tp = :C

        res(i+1) = ||F' * Z_i * Z_i' + Z_i * Z_i' * F + B' * B||F/||B' * B||F

    3. Note that all stopping criteria are checked only after a step
    with a real parameter or a "double step" with a pair of conjugate
    complex parameters. This ensures that Z*Z' is real, even if Z is
    not.

References:

    The algorithm is a slight modifivation of that proposed in:

    1. J.Li, F.Wang, and J.White.
    An efficient Lyapunov equation based approach for generating
    reduced-order models of interconnect.
    Proceedings of the 36th IEEE/ACM Design Automation Conference,
    New Orleans, LA, 1999.

    Another (though more expensive) low rank algorithm is proposed in:

    2. T. Penzl.
    A cyclic low rank Smith method for large sparse Lyapunov equations.
    To appear in SIAM Sci. Comp.

    See also:

    3. P. Benner, J. Li, and T. Penzl
    Numerical solution of large Lyapunov equations, Riccati equations,
    and linear-quadratic optimal control problems.
    In preparation.

    4. T. Penzl.
    LYAPACK (Users' Guide - Version 1.0).
    1999.

Internal remarks:

    Input data not completely checked!

    The procedure to generate real factors in case of complex parameters is
    different from that in [3]!

    The matrices SMi (i = 1:length(p)) for the "Sherman-Morrison trick"
    (only used if Bf and Kf are nonzero) are computed a priori, which
    is good in view of computation if parameters p(i) are used cyclically,
    but may be sometimes bad in view of memory demand, in particular, when
    length(p) is large.

    The stopping criterion related to the input parameter with_rs
    corresponds to the stagnation of the residual curve caused by
    round-off errors. Its performance depends on the constants stcf and
    min_rs. The iteration is stopped, as soon as

    (ra-rb) * (i-stcf+1) / ((r(1)-ra) * stcf) < min_rs

    and r(1)-ra > 0 and i >= stcf hold for stcf consecutive iteration
    steps. Here res(i+1) is the residual norm after the i-th LRCF-ADI step,
    r(i+1) = log(res(i+1)), ra = min(r(1:i-stcf+1)), rb = min(r(i-stcf+2:i+1)).

    stcf is also the number of consecutive steps, for which the criterion
    w.r.t. min_in must be fulfilled.
"""
function lp_lradi(A,B;p=[],Bf=[],Kf=[],K=[],max_it=100,tp=:B,rc=:C,
  min_res=1e-6,with_rs=:S,min_in=0,infoV=true)

tic()
isempty(p) && (p=lp_para(A))
const stcf = 10
const min_rs = .1

(tp!=:B && tp!=:C) && error("tp must be either :B or :C.")
(rc!=:R && rc!=:C) && error("rc must be either :R or :C.")

const tpB = (tp == :B)

const with_min_rs = (with_rs==:S)
const with_norm = (min_res>0)||with_min_rs
const make_real = (rc==:R)

const with_min_in = min_in>0

const with_BK = !isempty(Bf)

const l = length(p)

if tpB
  const n,m = size(B)
else
  const m,n = size(B)
end

const Ide = speye(size(A,1))
LP_L = Array(Complex{Float64},size(A)...,l)
LP_U = Array(Complex{Float64},size(A)...,l)
for i = 1:l
  LP_L[:,:,i],LP_U[:,:,i],~=lu(full(A+p[i]*Ide),Val{false})
end

if with_BK
  SM=Array(Complex{Float64},size(Bf)...,l)
  Im = speye(size(Bf,2))
  if tpB
    # SMi = TM*inv(I-Kf'*TM),
    # where TM = inv(F+p(i)*I)*Bf
    # (These are the columns of the LR terms for the
    # rank corrections of the "inverse"
    # in the Sherman-Morrison formulae.)

    for i = 1:l
      TM = LP_U[:,:,i]\(LP_L[:,:,i]\Bf)
      SM[:,:,i] = TM/(Im-Kf'*TM)
    end
  else  # tp==:C
    # SMi = TM*inv(I-Bf'*TM),
    # where TM = inv(F.'+p(i)*I)*Kf
    # (These are the columns of the LR terms for the
    # rank corrections of the "inverse"
    # in the Sherman-Morrison formulae.)
    for i = 1:l
      TM = LP_L[:,:,i].'\(LP_U[:,:,i].'\Kf)
      SM[:,:,i] = TM/(Im-Bf'*TM)
    end
  end
end

# Initialize QR factorization for norm computation
if with_norm
  res0,nrmQ,nrmR,nrmbs = lp_nrmu(A,B,Bf,Kf,tpB,[],[],[],[])
  res = [1.]
  res_log = [log(res0)]
  # Vector containing log of residual norms
  # corresponds to r(:) in prolog.
end

flag = "I"
if with_min_in
  nrmF_Z_2 = 0
  nrmF_rec = +Inf*ones(stcf,1)
  # Current squared Frobenius norm of Z
  # Records the values of
  # ||Z_nc||_F^2 / ||Z_new||_F^2 (see prolog)
  # for the last stcf iteration steps.
end


i_p = 1;                   # Pointer to i-th entry in p(:)
is_compl = imag(p[1])!=0;  # is_compl = (current parameter is complex.)
is_first = true;           # is_first = (current parameter is the first
                           #            of a pair, PROVIDED THAT is_compl.)


local V, Z
for i = 1:max_it       # The iteration itself

  if i==1

    if tpB
      # V = last columns of Cholesky factor Z
      # Initialize:
      # V := sqrt(-2*real(p[1]))*inv(F+p[1]*I)*B

      V=LP_U[:,:,1]\(LP_L[:,:,1]\B)
      with_BK && (V = V+SM[:,:,1]*(Kf'*V))
      V = sqrt(-2*real(p[1]))*V

    else #( tp = :C )
      # Initialize:
      # V := sqrt(-2*real(p[1]))*inv(F.'+p[1]*I)*B'

      V = LP_L[:,:,1].'\(LP_U[:,:,1].'\B')
      with_BK && (V = V+SM[:,:,1]*(Bf'*V))
      V = sqrt(-2*real(p[1]))*V
    end

    Z = V;             # Note:  Z*Z' = current ADI iterate

  else # i > 1

    p_old = p[i_p]

    i_p = i_p+1
    i_p>l && (i_p = 1)    # update current parameter index

    if is_compl && is_first
      is_first = false
      if i_p==1
        error("Parameter sequence ends in the 'middle' of a complex pair!")
      end
      if p[i_p]!=conj(p_old)
        error("Parameters p[i] must be either real or pairs of conjugate complex numbers.")
      end
    else
      is_compl = imag(p[i_p])!=0
      is_first = true
    end

    if tpB
      # Evaluate
      #   V := sqrt(real(p[i_p])/real(p_old))*...
      #        (V-(p[i_p]+conj(p_old))*inv(F+p[i_p]*I)*V)

      TM = LP_U[:,:,i_p]\(LP_L[:,:,i_p]\V)
      with_BK && (TM = TM+SM[:,:,i_p]*(Kf'*TM))
      TM = V-(p[i_p]+conj(p_old))*TM
      V = sqrt(real(p[i_p])/real(p_old))*TM

    else   #  tp==:C

      # Evaluate
      #   V := sqrt(real(p[i_p])/real(p_old))*...
      #        (V-(p[i_p]+conj(p_old))* inv(F.'+p[i_p]*I)*V)

      TM = LP_L[:,:,i_p].'\(LP_U[:,:,i_p].'\V)
      with_BK && (TM = TM+SM[:,:,i_p]*(Bf'*TM))
      TM = V-(p[i_p]+conj(p_old))*TM
      V = sqrt(real(p[i_p])/real(p_old))*TM
    end


    # Form new iterate Z.

    !is_compl && (V = real(V))

    Z = [Z V]

    # Make last 2*m columns real.
    if make_real && is_compl && !is_first
      for j = (i-1)*m+1:i*m
        U1,S1,V1 = svd([real(Z[:,j-m]) real(Z[:,j]) imag(Z[:,j-m]) imag(Z[:,j])])
        S1 = diagm(S1)
        U2,S2,V2 = svd(V1[1:2,1:2]'*V1[1:2,1:2]+V1[3:4,1:2]'*V1[3:4,1:2],thin=false)
        TMP = U1[:,1:2]*S1[1:2,1:2]*U2*diagm(sqrt(S2))
        Z[:,j-m] = TMP[:,1]
        Z[:,j] = TMP[:,2]
      end

    end

  end

  if with_norm
    # Compute residual norm

    resnrm,nrmQ,nrmR,nrmbs = lp_nrmu( A, B, Bf, Kf, tpB, V, nrmQ, nrmR, nrmbs)
    push!(res_log, log(resnrm))
    akt_res = resnrm/res0
    push!(res, akt_res)

    infoV && println("LRCF-ADI step $i -- norm. residual = $akt_res")

    # After pair of complex parameters or
    # a real parameter, check stopping criteria
    # based on residual norm.
    if  !(is_compl && is_first)

      if akt_res <= min_res
        flag = "R"
        break
      end

      if with_min_rs && i>=2*stcf
        ra = minimum(res_log[1:i-stcf+1])
        rb = minimum(res_log[i-stcf+2:i+1])
        if res_log[1]-ra > 0 && (ra-rb)*(i-stcf+1)/((res_log[1]-ra)*stcf) < min_rs
          flag = "S"
          break
        end
      end

    end

  end
  # Check stopping criteria based on increase in ||Z_i||_F.
  if with_min_in
    nrmF_V_2 = sum(sum(abs(V).*abs(V)));    # Note the "abs"; V is complex.
    nrmF_Z_2 = nrmF_Z_2 + nrmF_V_2
    nrmF_rec[1:stcf-1] = nrmF_rec[2:stcf]
    nrmF_rec[stcf] = nrmF_V_2/nrmF_Z_2
    if !(is_compl && is_first) && i>stcf && all(nrmF_rec .< min_in)
      flag = "N"
      break
    end
  end
end

infoV && println("Computational Time:  ", toq())

Z, flag, res
end

"""
Based on lp_nrmu.m from LYAPACK 1.0

Helper function for lp_lradi

Using updated QR factorizations, this routine computes efficiently
a sequence of norms which correspond to either of the following types:

for tp = :B

    nrm = || F*Z_i*Z_i' + Z_i*Z_i'* F' + B*B' ||F,

for tp = :C

    nrm = || F' * Z_i * Z_i' + Z_i * Z_i' * F + B' * B ||F.

Here, F = A-Bf*Kf'.

The matrices Z_i must have much more rows than columns and they
have to obey the recursion

    Z_i = [ Z_{i-1}  V ].

Calling sequence:

    nrm, nrmQ, nrmR, nrmbs = ...
    lp_nrmu(  A, B, Bf, Kf, tpB, V, nrmQ, nrmR, nrmbs )

Input:

    tpB        (= true or false) the type of the norm

    Bf        real matrix Bf
              Set Bf = [] if not existing or zero!

    Kf        real matrix Kf
              Set Kf = [] if not existing or zero!

    B         n-x-m or m-x-n matrix B (must be real)

    V         n-x-r matrix V (may be complex)

    nrmQ,

    nrmR,

    nrmbk     variables for internal use (they contain the data of
              the updated QR factorization).

Output:

    nrm       the current value of the Frobenius norm nrm_i.

    nrmQ,

    nrmR,

    nrmbk     variables for internal use (they must be output
              parameters).

Remarks:

    Using this routine for computing the norms within an iteration
    can save a lot of computation compared to repeatedly calling
    'lp_nrm'.

    This routine must be used as follows:

    1. Before the iteration starts, it must be invoked with the
       parameters V = nrmQ = nrmR = nrmbl = []. This initializes the
       variables nrmQ, nrmR, nrmbl.

    2. The routine must be called in each step of the iteration:

Example:

    ...

    nrm,nrmQ,nrmR,nrmbs = lp_nrmu(A,B,[],[],tpB,[],[],[],[])

    ...

    Z = zeros(n,0)

      ...

      for i = 1:100

        ...

        V = ...

        ...

        Z = [ Z V ];    # iteration:  Z_i = [ Z_i-1 V ]

        ...

        nrm,nrmQ,nrmR,nrmbs = lp_nrmu(A,B,[],[],tpB,V,nrmQ,nrmR,nrmbs)

        ...

      end
      ...
"""
function lp_nrmu( A, B, Bf, Kf, tpB, V, nrmQ, nrmR, nrmbs )

const with_BK = length(Bf)>0

const n= size(A,1);                    # Get system order.

if isempty(V)
  # The routine is called for the first time (before the iteration
  # 'outside' is started (i.e. V = zeros(n,0) or []!)) in order to
  # initialize the QR factorization nrmQ*nrmR correctly.

  if tpB
    nrmbs = [size(B,2)]
    nrmQ,nrmR = qr(B)
  else
    nrmbs = [size(B,1)]
    nrmQ,nrmR = qr(B')
  end
  nrm = vecnorm(nrmR*nrmR')

else
  # The routine is called during the iteration 'outside'.
  # The QR factorization is updated and the norm is computed.

  # Update of the QR factorization.
  push!(nrmbs,size(V,2))
  lw = size(nrmQ,2)
  lz = size(V,2)

  if tpB
    TM = A*V
    with_BK && (TM = TM-Bf*(Kf'*V))
    Z = [ TM V ]
  else
    TM = A'*V
    with_BK && (TM = TM-Kf*(Bf'*V))
    Z = [ TM V ]
  end

  for j = 1:2*lz
    a = Array(Complex{Float64},lw,1)
    #a = complex(zeros(lw,1))
    t = Z[:,j]
    for k = 1:lw
      u = nrmQ[:,k]
      alpha = dot(u,t)
      t = t-alpha*u
      a[k] = alpha[1]
    end
    beta = norm(t)
    nrmQ = [nrmQ t/beta]
    #nrmQ = push!(nrmQ, t/beta...)
    nrmR = [ nrmR a ; zeros(1,lw) beta ]
    lw = lw+1
  end

  # Computation of  nrmR * [permutation matrix] * nrmR'
  RT = copy(nrmR)
  ie2 = nrmbs[1]
  for j = 2:length(nrmbs)
    ia1 = ie2+1
    ie1 = ie2+nrmbs[j]
    ia2 = ie1+1
    ie2 = ie1+nrmbs[j]
    TMP = RT[:,ia1:ie1]
    RT[:,ia1:ie1] = RT[:,ia2:ie2]
    RT[:,ia2:ie2] = TMP
  end
  nrm = vecnorm(nrmR*RT')
end
nrm, nrmQ, nrmR, nrmbs
end


"""
    `lp_para(A,Bf,Kf,l0,kp,km,b0)`

Based on lp_para.m from LYAPACK 1.0

Estimation of suboptimal ADI shift parameters for the matrix

F = A - Bf * Kf'.

Calling sequence:

p = lp_para(A,Bf,Kf,l0,kp,km,b0)

Input:

    Bf        matrix Bf
              Set Bf = [] if not existing or zero!
    Kf        matrix Kf
              Set Kf = [] if not existing or zero!
    l0        desired number of shift parameters (kp+km > 2*l0)
              (The algorithm delivers either l0 or l0+1 parameters!)
    kp, km    numbers of Arnoldi steps w.r.t. F and inv(F),
              respectively (kp, km << order of A)
    b0         Arnoldi start vector (optional; chosen at random if not
              provided).

Output:

    p         an l0- or l0+1-vector of suboptimal ADI parameters

Remarks:

    Typical values are l0 = 10..40, kp = 20..80, km = 10..40.
    The harder the problem is the large values are necessary.
    Larger values mostly result in a faster convergence, but also in a
    larger memory requirement.
    However, for "well-conditioned" problems small values of l0 can
    lead to the optimal performance.

References:

    1. T. Penzl.
    LYAPACK (Users' Guide - Version 1.0).
    1999.

Input data not completely checked!
"""
function lp_para(A,Bf=[],Kf=[],l0=NaN,kp=NaN,km=NaN,b0=[])

err_code = 0
n=size(A,1); # Get system order.

if isnan(l0) && isnan(kp) && isnan(km)
  if n > 51
    kp=50;
    km=25;
    l0=20;
  elseif n > 15
    kp = n - 5
    km = Int(round(kp/2))
    l0 = km-2
  else
    kp = n-2
    km = Int(round(kp/2))
    l0 = km-1
  end
else
  kp >= n && error("kp must be smaller than n!")
  km >= n && error("km must be smaller than n!")
  2*l0 >= kp+km && error("2*l0 must be smaller than kp+km!")
end

isempty(b0) && (b0 = randn(n,1))
b0 = (1/norm(b0))*b0

rwp = []
rwm = []
rw = []
Hp = []
Hm = []

if kp > 0
  Hp,V = lp_arn(A,Bf,Kf,kp,b0,:p)
  rwp = eigvals(Hp[1:kp,1:kp]);                  # =: R_+
  rw = rwp
end

if km > 0
  Hm,V = lp_arn(A,Bf,Kf,km,b0,:m)
  rwm = ones(km,1)./eigvals(Hm[1:km,1:km]);      # =: 1 / R_-
  push!(rw,rwm...);                           # =: R
end

if any(real(rw) .>= zeros(size(rw)))
  err_code = 1
  println("These are the Ritz values computed by the Arnoldi process w.r.t. F: $rwp
  These are the Ritz values computed by the Arnoldi process w.r.t. inv(F): $rwm

  ####################################################################
  WARNING in 'lp_para': NON-STABLE RITZ VALUES DETECTED!!!

  This is quite a serious problem, that can be caused by
  (i)   non-stable matrices F (Be sure that F is stable. ADI like
        methods only work for stable or antistable problems. If your
        Lyapunov equation is antistable, multiply it by -1.)
  (ii)  matrices F that are stable but have an indefinite symmetric
        part (This is THE weak point of this algorithm. Try to work
        with the 'reduced' Ritz values, i.e., the unstable values are
        simply removed. This is not an elegant measure but it may work.
        However, the convergence of ADI can be poor. This measure is
        taken automatically. Another measure might be to enlarge the
        values of kp or km, and run the program again.
  (iii) matrices F with a negative definite, but ill-conditioned
        symmetric part (This is quite unlikely. The problem is
        caused by round-off errors).

  #####################################################################


  NOTE: The unstable Ritz values will be ignored in the further computation!!! ")
  rw0 = copy(rw)
  rw = []
  for j = 1:length(rw0)
    if real(rw0[j])<0
      push!(rw,rw0[j])
    end
  end
end

p = lp_mnmx(rw,l0)

p
end


"""
Based on lp_arn_m.m and lp_arn_p.m from LYAPACK 1.0 and includes many of the
same comments much of the same description.

helper function for lp_para

Arnoldi method w.r.t. inv(F), where F = A-Bf*Kf'.

Calling sequence:

    H,V = lp_arn_m(A,Bf,Kf,k,r)

Input:

    A         Matrix A

    Bf        matrix Bf
              Set Bf = [] if not existing or zero!

    Kf        matrix Kf

              Set Kf = [] if not existing or zero!

    k         number of Arnoldi steps (usually k << n, where
              n is the order of the system)

    r         initial n-vector
              (optional - chosen by random, if []).

    pm        %TODO write description

Output:

    H         matrix H ((k+1)-x-k matrix, upper Hessenberg)

    V         matrix V (n-x-(k+1) matrix, orthogonal columns).


Method:

    The ("inverse") Arnoldi method produces matrices V and H such that

    V(:,1) in span{r},
    V' * V = eye(k+1),
      inv(F) * V(:,1:k) = V * H.

Remark:

    This implementation does not check for (near-)breakdown!

"""
function lp_arn(A,Bf,Kf,k,r,pm)


const with_BK = !isempty(Bf)

const n = size(A,1)                 # Get system order.
k >= n-1 && error("k must be smaller than the order of A!")
isempty(r) && (r = randn(n,1))

V = zeros(n,k+1)
H = zeros(k+1,k)

V[:,1] = (1.0/norm(r))*r

beta = 0

LP_L,LP_U,~ = lu(full(A),Val{false})

if with_BK && pm==:m
  # SM = inv(F)*Bf*inv(I-Kf'*inv(F)*Bf)
  # (This is the main part of the term needed for the low
  # rank correction in the Sherman-Morrison formula.)
  Im = speye(size(Bf,2))
  TM = LP_U\(LP_L\Bf)
  SM = TM/(Im-Kf'*TM)
end

for j = 1:k

  if j > 1
    H[j,j-1] = beta
    V[:,j] = (1.0/beta)*r
  end

  if (pm == :m)
    w = LP_U\(LP_L\V[:,j])
    with_BK && (w = w + SM*(Kf'*w)); # LR correction by SM formula
  else
    w = A*V[:,j]
    with_BK && (w = w-Bf*(Kf'*V[:,j]))
  end
  r = w

  for i = 1:j
    H[i,j] = dot(V[:,i],w)
    r = r-H[i,j]*V[:,i]
  end

  beta = norm(r)
  H[j+1,j] = beta

end

V[:,k+1] = (1.0/beta)*r

H,V
end


"""
Based on lp_mnmx.m from LYAPACK 1.0 and includes many of the same comments and
much of the same description.

Suboptimal solution of the ADI minimax problem. The delivered parameter
set is closed under complex conjugation.

Calling sequence:

    p = lp_mnmx(rw,l0)

Input:

    rw        a vector containing numbers in the open left half plane, which
              approximate the spectrum of the corresponding matrix, e.g.,
              a set of Ritz values. The set must be closed w.r.t. complex
              conjugation

    l0        desired number of shift parameters (length(rw) >= l0)
              (The algorithm delivers either l0 or l0+1 parameters!).

Output:

    p         an l0- or l0+1-vector of suboptimal ADI parameters

"""
function lp_mnmx(rw,l0)

if length(rw)<l0
  error("length(rw) must be at least l0.")
end

max_rr = +Inf;                       # Choose initial parameter (pair)
p0 = 0
for i = 1:length(rw)
  max_r,~ = lp_s(rw[i],rw)
  if max_r < max_rr
    p0 = rw[i]
    max_rr = max_r
  end
end

imag(p0)!=0 ? (p = [ p0; conj(p0) ]) : (p = [p0])

max_r,i = lp_s(p,rw);         # Choose further parameters.

while size(p,1) < l0

  p0 = rw[i]
  imag(p0)!=0 ? push!(p,p0,conj(p0)) : push!(p,p0)

  max_r,i = lp_s(p,rw)

end

p
end


"""
Based on lp_s.m from LYAPACK 1.0 and includes much of the same description

Helper function for lp_para

Computation of the maximal magnitude of the rational ADI function over
a discrete subset of the left complex half plane.

Calling sequence:

    max_r,ind = lp_s(p,set)

Input:

    p        vector of ADI parameters

    set      vector representing the discrete set.

Output:

    max_r    maximal magnitude of the rational ADI function over set

    ind      index - maximum is attained for set(ind).
"""
function lp_s(p,set)

max_r = -1
ind = 0

for i = 1:length(set)

  x = set[i]

  rr = 1
  for j = 1:length(p)

    rr = rr*abs(p[j]-x)/abs(p[j]+x)

  end

  if rr > max_r

    max_r = rr
    ind = i

  end

end
max_r,ind
end


# end module matrixEqs
end
