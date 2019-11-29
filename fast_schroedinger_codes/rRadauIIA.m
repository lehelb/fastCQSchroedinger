function r = rRadauIIA(z,RK) 

% 1 + z*b^T*(Id - z*A)^(-1)*e  
if (RK == 1)
    r = 1+z./(1-z); 
elseif (RK == 2)
    r = 2.*(3+z)./(6-4.*z+z.^2); 
elseif (RK == 3)
    r = -3.*(20+8.*z+z.^2)./(-60+36.*z-9.*z.^2+z.^3); 
else
    error('Only 1,2,and 3 stage Radau IIA methods are implemented')
end
