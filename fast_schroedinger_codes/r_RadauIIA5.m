function r = r_RadauIIA5(z) 

% 1 + z*b^T*(Id - z*A)^(-1)*e  

r = -3.*(20+8.*z+z.^2)./(-60+36.*z-9.*z.^2+z.^3); 
