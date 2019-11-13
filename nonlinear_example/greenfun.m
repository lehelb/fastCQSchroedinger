function out=greenfun(lambda,x)
%lam:scalar
%x:array

out=exp(-sqrt(lambda)*abs(x))./(2*sqrt(lambda));
