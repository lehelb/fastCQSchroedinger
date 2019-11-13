function f = heaviside(t)

f = ones(size(t));
I = t < 0;
f(I) = 0;