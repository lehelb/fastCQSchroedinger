function Delta = DeltaRadauIIA3(z) 

% Delta(z) = (A + z/(1-z)*e*b^T)^-1  
% thank to maple
Delta =  [3./2,1./2-2.*z;
         -9./2,5./2+2.*z]; 
