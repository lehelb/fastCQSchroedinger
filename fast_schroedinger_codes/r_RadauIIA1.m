function r = r_RadauIIA1(z) 

% 1 + z*b^T*(Id - z*A)^(-1)*e  

r = 1+z./(1-z); 
