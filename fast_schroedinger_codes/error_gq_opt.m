function f = error_gq_opt(Q,B,Ljm1,d,RK,dt,n0)

if (RK == 1) %BDF1     
    rn = @(z) 1./abs(1-z); qn = @(z) 1./abs(1-z);
else %currently just 2-stage Radau       
    rn = @(z) abs((2*z+6)./(z.^2-4*z+6)); qn = @(z) sqrt((9/2)^2+abs(3/2-z).^2)./abs(z.^2-4*z+6); 
    if (RK > 2)
        disp(['Choice of parameters for 2-stage Radau IIA used. Might not be ' ...
              'optimal for the 3-stage version']);
    end
end

M1 = 20; M2 = 20;
rho_max = 1+(2/B)*(1+sqrt(1+B));
rhos = linspace(1,rho_max,M1);
rhos = rhos(2:end);

f = inf;
th = linspace(0,pi,M2);
for j = 1:length(rhos)
    rho = rhos(j);
    z = 1+(B/4)*(rho*exp(1i*th)+(1/rho)*exp(-1i*th)+1);    
    hnj = (1./sqrt(abs(z))).*exp(d*sqrt(Ljm1)*sqrt(abs(z))).*...
        abs(qn(-dt*z).*(rn(-dt*z).^(n0+1)));            
    f = min(f,(2*B*sqrt(Ljm1)*dt/pi)*(rho^(-2*Q+1)/(rho-1))*max(hnj));
end