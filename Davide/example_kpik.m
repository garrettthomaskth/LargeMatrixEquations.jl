%close all
%clear all
%clc

nh=3;
n=nh^2;
T=diag(2*ones(nh,1))+diag(-ones(nh-1,1),1)+diag(-ones(nh-1,1),-1);
I=speye(nh);
A=-(kron(T,I)+kron(I,T));
%E=speye(n);
%LE=E;
B=randn(n,1);
%B=0.2:0.1:1;
%B=B';
m=100;
tol=1e-9;
tolY=1e-12;
[Z,r]=kpik(A,B,m,tol,tolY);

fprintf('final true relative residual norm: \n')
disp(norm((A*Z)*Z'+Z*(Z'*A')+B*B','fro')/norm(B,'fro')^2)    %this matrix should never be formed for n large 
