function [X,W,n0,Q1,Q2s]=quadrature_cqw(d,tol,T,dt,B,RK)
Qmax = 200;
% weights 0, \dots, n_0 computed directly and n_0+1 etc use quadrature


%constants for estimates
if (RK == 1) %BDF1     
    b = .5; nu = -log(1-b)/b; Cq = exp(nu*b); CA = 1; %From the paper ffde
    rn = @(z) 1./abs(1-z); qn = @(z) 1./abs(1-z);
else %currently just 2-stage Radau    
    b = 1; nu = 1.0735; Cq = 1.4733; %From the paper ffde
    CA = 2.1213;%need to check what Cq1 is for Radau
    rn = @(z) abs((2*z+6)./(z.^2-4*z+6)); qn = @(z) sqrt((9/2)^2+abs(3/2-z).^2)./abs(z.^2-4*z+6); 
end




A = 1;%??
xi=A/sqrt(dt); %Not so big. Requires bigger n0
%xi=A/dt;
%xi = A*dt^(-3/4);
%xi = 2;
n0 = 1;
err = inf;
while (err > tol/2)
    n0 = n0+1;
    %err = cosh(d*sqrt(xi/2))*(dt/(2*pi))*integral(@(y) abs(qn(-dt*xi+1i*y*dt).*(rn(-xi*dt+1i*y*dt).^(n0+1)))...
    %        .*(xi^2+y.^2).^(-.25),0,inf);        
    err1 = (dt/(4*pi))*integral(@(y) exp(-d*real(exp(-1i*pi/4)*sqrt(-xi+1i*y))).*abs(qn(-dt*xi+1i*y*dt).*(rn(-xi*dt+1i*y*dt).^(n0+1)))...
            .*(xi^2+y.^2).^(-.25),0,inf);        
    err2 = (dt/(4*pi))*integral(@(y) exp(-d*real(exp(-1i*pi/4)*sqrt(-xi-1i*y))).*abs(qn(-dt*xi+1i*y*dt).*(rn(-xi*dt+1i*y*dt).^(n0+1)))...
            .*(xi^2+y.^2).^(-.25),0,inf);            
    err = err1+err2;
end
n0_est=ceil((d*sqrt(xi/2)-log(4*xi*(gamma(3/4)^2)*tol/(CA*sqrt(2*pi))))/(test_gamma(xi*dt,RK)*xi*dt));  %to achieve a truncation error bounded by tol. Follows from Prop 9.
[n0 n0_est]
% Do Gauss-Jacobi on interval [0,L0]
%L0 = min(1/T,1/d^2); %L0=1/T %As in ffde
L0 = 2/T;
%number of subintervals
J = round(log(xi/L0)/log(1+B));
B = (xi/L0)^(1/J)-1;

%as there are J intervals
tol = (tol/2)/J;


%G-J interval
rho_max = 1+2*b/(L0*dt)+sqrt((2*b/(L0*dt))^2+4*b/(L0*dt));
err = inf;
Q1 = 0;


while (Q1 < Qmax && err > tol)
    Q1 = Q1+1;
    rho_opt = 2*Q1/(d*sqrt(3*L0/2)+nu*T*L0/2) + sqrt(1+(2*Q1/(d*sqrt(3*L0/2)+nu*T*L0/2))^2);
    if (rho_opt < rho_max) && (rho_opt > 2+sqrt(3))
        err = Cq*(4*dt/pi)*sqrt(L0)*(rho_opt/(rho_opt-1))*(exp(1)*(d*sqrt(3*L0/2)+nu*T*L0/2)/(4*Q1))^(2*Q1);        
    else
        err = Cq*(4*dt/pi)*sqrt(L0)*exp((d*sqrt(6/L0)+nu*T)*b/dt)/(rho_max-1)*rho_max^(-2*Q1+1);                
    end
end

if (rho_opt > rho_max) || (rho_opt <= 2+sqrt(3))
    disp('rho_max used in GJ')
    disp([rho_opt rho_max])
end

Q2s = zeros(J,1);
%Computation of all quadrature weights and nodes 

%Nodes and weights for the integral in [0,L0]
[xjac,wjac] = gaussj(Q1,0,-1/2); 
wjac = wjac.*sqrt(1+xjac)*L0/(4*pi*1i); %Includes correction in the integrand 
xjac = (xjac+1)*L0/2; 
X=xjac; W=wjac;

%Nodes and weights in [L_{j-1},L_j], j=1,...,J
h = waitbar(0,'Computing Gauss intervals');
for j = 0:(J-1)
    Lj = L0*((1+B)^j); %L_{j-1}
    if j==J-1 %avoids going beyond xi with the last subinterval
        B=xi/Lj-1;
    end
   
    Q2 = 0;
    err = inf;
    while (Q2 < Qmax && err > tol)
        Q2 = Q2+1;
        err = error_gq_opt(Q2,B,Lj,d,RK,dt,n0);
    end    
    Q2s(j+1) = Q2;
    [xg,wg] = gauss(Q2);
    
    xgj = (xg+1)*B*Lj/2+Lj;
    wgj = wg.'*B*Lj/(4*pi*1i);
    X=[X;xgj];
    W=[W;wgj];
    waitbar((j+1)/J,h)
end
close(h)
