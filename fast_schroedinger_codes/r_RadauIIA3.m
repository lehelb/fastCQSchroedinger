function r = r_RadauIIA3(z) 

% 1 + z*b^T*(Id - z*A)^(-1)*e  

r = 2.*(3+z)./(6-4.*z+z.^2); 
