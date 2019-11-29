function [A,b,c] = RKdata(RK)
%function [A,b,c] = RKdata(RK)
%returns coefficients in Butcher notation
%RK = 0 for backward Euler
%RK = 1 for 2 stage Radau IIA
%RK = 2 for 3 stage Radau IIA

if (RK == 1)
    %BDF1
    A = 1; c = 1; b = 1;        
elseif (RK==2)
    %radau IIa of order 3 stage order 2
    A = [5/12 -1/12; 3/4 1/4]; c = [1/3; 1]; 
    b = [3/4; 1/4];        
elseif (RK==3)
%radau IIa order 5, stage order 3
    A = [(88-7*sqrt(6))/360 (296-169*sqrt(6))/1800 (-2+3*sqrt(6))/225;
         (296+169*sqrt(6))/1800 (88+7*sqrt(6))/360 (-2-3*sqrt(6))/225;
        (16-sqrt(6))/36 (16+sqrt(6))/36 1/9];
    b = [(16-sqrt(6))/36 (16+sqrt(6))/36 1/9].';
    c = [(4-sqrt(6))/10 ; (4+sqrt(6))/10 ; 1];        
else
    error('only RK types 1,2,3 possible, see help\n');    
end
