function [x,w] = gaussjb(n,alp)
%function [x,w] = gaussjb(n,alp)
% Nodes and weights for Gauss integration on interval [-1,1]
% with weight (1+x)^alp   
% Borrows from gaussj.m in the Matlab SC Tooolbox and in turn from GAUSSJ 
% in the SCPACK Fortran

a(1) = alp/(alp+2);
b(1) = sqrt(4*(1+alp)/((alp+3)*(alp+2)^2));
N = 2:n;
a(N) = alp^2 ./ ((alp+2*N).*(alp+2*N-2));
N = 2:(n-1);
b(N) = sqrt(4*N.^2.*(N+alp).^2./(((alp+2*N).^2-1).*(alp+2*N).^2));
if n > 1
  [V,D] = eig(diag(a) + diag(b,1) + diag(b,-1));
else
  V = 1; D = a;
end
x = diag(D);w = (2^(alp+1)/(alp+1))*(V(1,:)').^2;
[x,ind] = sort(x); w = w(ind);