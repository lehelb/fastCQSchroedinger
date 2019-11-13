function [solminusa,solplusa]=freesol(t,a,lf,le,Nf,Ne,alfa,beta)
%Solution to the free Schr�dinger equation for the initial data:
% psi_0(x)=alfa*phi_f(x)+beta*phi_e(x), x in R.

IA = @(la,t) pi/sqrt(la)*exp(1i*la*t).*(1-erfz(sqrt(1i*la*t)));

IB = @(la,t) pi/(2*sqrt(la))*exp(1i*la*t).*(exp(2*a*sqrt(la))*(1-erfz(sqrt(1i*la*t)+a./sqrt(1i*t)))...
    +exp(-2*a*sqrt(la))*(1-erfz(sqrt(1i*la*t)-a./sqrt(1i*t))));

IAlet = IA(le,t);IAlft = IA(lf,t);
IBlet = IB(le,t);IBlft = IB(lf,t);

solminusa=alfa*Nf*(IAlft+IBlft)+beta*Ne*(IAlet-IBlet);
solminusa=solminusa/(2*pi);

solplusa=alfa*Nf*(IAlft+IBlft)-beta*Ne*(IAlet-IBlet);
solplusa=solplusa/(2*pi);