function rksm(A,E,EL,B,m,tol,s1,emax,ch,tolY)
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
# ch      ch=1  complex poles  ch=0 real poles
# tolY    truncation tolerance for final solution, e.g., tolY=1e-12
#
# Output:
#
# Z    factor of approximate solution  X = Z Z'
#      If ch=1,  Z may be complex
# nrmrestot  history of residuals
#
# Hints:
# 1) Before the call permute entries of A, E and B so as to
#    limit fill-in in the system solves
# 2) Provide "comfortable" (loose bounds) estimates s1, emax
#
#
#  Please contact V. Simoncini for any problem you may encouter when
#  running the code
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


n=size(A,1)
B=full(B)
p=size(B,2)
I=speye(p)
O=0*I
uno=ones(1,p)
Lres=EL\B
V,irr=qr(Lres)
rr=inv(irr)
nrmb=vecnorm(inv(rr))^2
beta=V'*Lres
beta2=beta*beta'
s=Complex{Float64}[]
print("     no its     backward error\n")
VV=zeros(n,p*(m+2))
VV[1:n,1:p]=V
H=zeros(p*(m+2),p*(m+1))
nrmrestot=[]
nrma=vecnorm(A)

complexType = Union{Array{Complex{Int},2},Array{Complex{Bool},2},Array{Complex{Float64},2}}
isCmplx = false
if (isa(A,complexType))
  H = complex(H)
  VV = complex(VV)
  beta2 = complex(beta2)
  isCmplx = true
end

if (vecnorm(E-speye(n))>1e-14)
  condestE=cond(E);
  singE=condestE/vecnorm(E);
else
  singE=1;
end

if (norm(A-E-(A-E)',1)<1e-14)
  symm=true
else
  symm=false
end

newAv=EL\(A*(EL'\V));
K=V'*newAv;
push!(s,s1); #s[1]=s1;
eH=eig(K)
eHpoints = sort([s1,emax])
snew=newpolei(eHpoints,eH[1],s1*uno');

if real(snew)<0
   snew=-real(snew)+im*imag(snew);
end
push!(s,snew); #s[2]=snew;

# additional steps
cmplxflag=false;

i=0;
global i1, j1, js, j1s, Y
while i < m

  i=i+1;

  paired=0;
  while (paired==0)

    i1=i+1;
    w=EL*V;
    wrk = (A-snew*E)\w;
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
  d=-VV[1:n,1:js]*(Y*(H[1:ih*p,1:ih*p]'\[zeros(p*(ih-1),p);I])*H[p*ih+1:p*ih1,p*ih-p+1:p*ih]');
  ### Check
  U=[-V*s[end]  d u1 ];
  extra,rr=qr(full(U));
  nrmres=vecnorm(rr*sparse([O I O; I O I; O I O ])*rr')/(nrmb+singE*nrma*nrmx);
  nrmrestot=[nrmrestot; nrmres];

  println([i,nrmres])

  if (nrmres<tol)
    break
  end


  # New poles and zeros
  eH=eig(K)[1];
  eH=eH[sortpermComplex(eH)]

  # may cause problems
  eHorig=eH+0;

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
      else                                  # if all real eigs, no convex hull possible
        eHpoints = sort([s1; emax.';-real(eH)]);
      end


   else   # Real poles s from real set. Compute complex roots of r_m via Ritz convex hull
       if (countnz(imag(eH))>0 && length(eH)>2)    # Roots lambdas come from convex hull too
         eH=full([eH;-s1;-emax.']);
         #ij=convhull(real(eH),imag(eH));
         #print("CONVEXHULL")
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


   gs=kron(s[2:i+1].',complex(uno))';


   snew = newpolei(eHpoints,eH,gs);
   if real(snew)<0
     snew=-real(snew)+im*imag(snew);
   end  #safeguard strategy


   # If pole is complex, include its conjugate
   if (imag(snew) !=0)
     cmplxflag=true;
   end
   #s[i+2]=snew;

   push!(s,snew)

   g1 = g;
   g2 = V'*(EL\(A*(EL'\VV[1:n,1:js])));
   g3 = V'*(EL\(A*(EL'\V)));
   ### Check
   K = [K g1; g2 g3];
   VV[1:n,js+1:j1s]=V;

end;

# Done
# Reduce rank of solution, if needed
sY,uY=eig(Y)
id=sortpermComplex(sY)
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
function ratfun(x,eH,s)
r = complex(zeros(length(x),1))
for j=1:length(x)
  xj = x[j]*complex(ones(length(s)))
  r[j]=abs(prod( (xj-s)./(xj-eH) ));
end
return r
end

##################################
function newpolei(eHpoints,eH,s)
snew=complex(zeros(length(eHpoints)-1,1))
for j=1:length(eHpoints)-1

    sval=linspaceComplex(eHpoints[j],eHpoints[j+1],200);

    sf = maximum(abs(ratfun(sval,eH,s)));
    jx = indmax(abs(ratfun(sval,eH,s)))

    snew[j]=sval[jx];
end
sn=maximum(abs(ratfun(snew,eH,s)));
jx=indmax(abs(ratfun(snew,eH,s)))
snew=snew[jx];
return snew
end

####################################
function sortpermComplex(arr)
  if (countnz(imag(arr)) == 0)
    return sortperm(real(arr))
  else
    return sortperm(abs(arr))
  end
end

##########################################
function linspaceComplex(start,stop,steps)
  if imag(start) == 0 && imag(stop) == 0
    return(linspace(start,stop,steps))
  else
    realSp = linspace(real(start),real(stop),steps)
    imagSp = linspace(imag(start),imag(stop),steps)
    return(realSp+im*imagSp)
  end
end

#######################################
function convhull(pnts)
    # Function to compute 2-D convex hull
    T = eltype(pnts) #get point type
    N = length(pnts) #number of pnts
    # sort the points lexicographically.
    # copy points into mutable container
    #pnts  = T[]
    #for x in points
    #    push!(pnts, x)
    #end

    #trivial case 0 or 1 point
    length(pnts) <= 1 && (return pnts)

    #=
    sort points lexicographically
    =#
    sort!(pnts, lt=lt2d)

    #=
    function to orient boundry using 2D cross product of OA and OB vectors,
     i.e. z-component of their 3D cross product.
     Returns a positive value, if OAB makes a counter-clockwise turn,
     negative for clockwise turn, and zero if the points are collinear.
    =#
    orient(b, pnt) = (real(b[end]-b[end-1]) * imag(pnt - b[end-1])) -
                      (imag(b[end-1] - b[end-1]) * real(pnt - b[end-1]))

    #build lower hull
    lower = T[]
    #buildhull!(lower, pnts, 1:1:N)
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

#=
less than comparator
=#
function lt2d(a, b)
    dx, dy =  real(a - b), imag(a - b);
    dx != 0 &&  (return <(dx, 0))
    return <(dy, 0)
end
