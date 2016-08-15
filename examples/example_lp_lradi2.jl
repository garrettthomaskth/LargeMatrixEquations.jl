#
#  This example is based on demo_m1.m from LAYPACK 1.0 and includes
#  many of the same comments.
#
#  MODEL REDUCTION BY THE ALGORITHMS LRSRM AND DSPMR. THE GOAL IS TO
#  GENERATE A REDUCED SYSTEM OF VERY SMALL ORDER.
#
#  This demo program shows the use of 'lp_lradi' and 'lp_para'

# -----------------------------------------------------------------------
# Generate test problem
# -----------------------------------------------------------------------
#
# This is an artificial test problem of a system, whose Bode plot shows
# "spires".
workspace()
include("../src/matrixEqs.jl")
using matrixEqs


A = speye(408,408); B = ones(408,1); C = ones(408,1);
A[1:2,1:2] = [-.01 -200; 200 .001];
A[3:4,3:4] = [-.2 -300; 300 -.1];
A[5:6,5:6] = [-.02 -500; 500 0];
A[7:8,7:8] = [-.01 -520; 520 -.01];
A[9:408,9:408] = diagm(-(1:400));

println("Problem dimensions:")

n = size(A,1);   # problem order (number of states)
m = size(B,2);  # number of inputs
q = size(C,1);   # number of outputs

# -----------------------------------------------------------------------
# Initialization/generation of data structures used in user-supplied
# functions and computation of ADI shift parameters
# -----------------------------------------------------------------------
#
# Note that A is a tridiagonal matrix. No preprocessing needs to be done.


println("Parameters for heuristic algorithm which computes ADI parameters:")
l0 = 10;   # desired number of distinct shift parameters
kp = 30;   # number of steps of Arnoldi process w.r.t. A
km = 15;   # number of steps of Arnoldi process w.r.t. inv(A)

b0 = ones(n,1);   # This is just one way to choose the Arnoldi start
                  # vector.
Bf = []#zeros(408,2)#
Kf = [];
# computation of ADI shift parameters
p = lp_para(A,Bf,Kf,l0,kp,km,b0);

println("Actual number of ADI shift parameters:");
l = length(p);

println("ADI shift parameters:");
p;


# -----------------------------------------------------------------------
# Solution of Lyapunov equations A*X+X*A' = -B*B' and
# A'*X+X*A = -C'*C
# -----------------------------------------------------------------------

println("Parameters for stopping criteria in LRCF-ADI iteration:")
max_it = 40;   # (will stop the iteration)
min_res = 1e-100;   # (avoided, but the residual history is shown)
with_rs = "N";   # (avoided)
min_in = 0;   # (avoided)

zk = "Z";
rc = "R";
Bf = ones(408,2)#[];
Kf = zeros(408,2)#[];
K = ones(408,408);
info = 2;

println("... solving A*XB+XB*A'' = - B*B''...");
tp = "B";

#figure(1), hold off; clf;   # (lp_lradi will plot residual history.)
#using ProfileView
ZB,flag,res =lp_lradi(A,B,p,max_it=200)
#Profile.clear()
#@profile lp_lradi(A,B,p,max_it=200)
#ProfileView.view()
#ZB,flag,res = lp_lradi(A,B,p)#,max_it,tp,zk,rc,min_res,with_rs,min_in,info);
# compute ZB

X=ZB*ZB'
println(norm(A*X+X*A'+B*B'))
#=
title("LRCF-ADI for CALE  AX_{B}+X_{B}A^T = -BB^T")
println("Termination flag:")
flag_B
println("Size of ZB:");
size_ZB = size(ZB)

println("... solving A''*XC+XC*A = - C''*C...");
tp = "C";

figure(2), hold off; clf;   # (lp_lradi will plot residual history.)

[ZC,flag_C] = lp_lradi(tp,zk,rc,name,Bf,Kf,C,p,max_it,min_res,...
              with_rs,min_in,info);
                                       # compute ZC

title("LRCF-ADI for CALE  A^T X_{C} + X_{C} A_ = -C^TC")
println("Termination flag:")
flag_C
println("Size of ZC:");
size_ZC = size(ZC)
=#
