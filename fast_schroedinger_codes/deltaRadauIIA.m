function delta = deltaRadauIIA(z,RK) 

% delta(z) = (A + z/(1-z)*e*b^T)^-1  
if (RK == 1)
    delta = 1-z;
elseif (RK == 2)
   delta =  [3./2,1./2-2.*z;
         -9./2,5./2+2.*z];  
elseif (RK == 3)
    nz = zeros(size(z));
    delta = cat(1,cat(1, cat(2,cat(2, nz+2+1./2.*6.^(1./2),nz-6./5+29./30.*6.^(1./2))   ,-6./5.*z-6./5.*z.*6.^(1./2)+2./5-4./15.*6.^(1./2)),...
                     cat(2,cat(2, nz-6./5-29./30.*6.^(1./2),nz+2-1./2.*6.^(1./2)),-6./5.*z+6./5.*z.*6.^(1./2)+2./5+4./15.*6.^(1./2))),...
               cat(2,cat(2,nz-1+8./3.*6.^(1./2),nz-1-8./3.*6.^(1./2)),-3.*z+5)); 
else
    error('Only 1,2,and 3 stage Radau IIA methods are implemented')
end
