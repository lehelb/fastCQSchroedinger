function [lf,le,Nf,Ne]=cnorm(gamma1,gamma2,a)
%Compute eigenvalues
%Good starting values for gamma1=gamma2=-.5, a=3
determinant = @(l) (1/gamma1+1./(2*sqrt(l))).*(1/gamma2+1./(2*sqrt(l)))...
        -exp(-4*sqrt(l)*a)./(4*l);

lf=fzero(determinant,.1);
le=fzero(determinant,.05);
 %Compute L2-norm of eigenfunctions
dx=1/100; x=-200:dx:200;
y = phifnn(a,lf,x).^2; y(1,end) = y(1,end)/2;
Nf = 1/sqrt(dx*sum(y));

z = phienn(a,le,x).^2; z(1,end)=z(1,end)/2;
Ne=1/sqrt(dx*sum(z));