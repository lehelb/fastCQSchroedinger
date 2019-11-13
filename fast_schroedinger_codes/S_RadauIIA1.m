function S = S_RadauIIA1(z) 

% S(z) = r^-1 * (I - z*A) *e* bT * (I - z*A)   

S = -1./(-1+z); 
