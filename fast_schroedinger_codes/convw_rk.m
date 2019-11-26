function [omega,omega0] = convw_rk(N,dt,RK,K)
%function f = convw_rk(N,dt,RK,Ks)
%returns conv. weight \omega_j(K) for j = 0, .., N
%RK = (1 RadauIIA1) (2 RadauIIA2) (3 RadauIIA3) (4 LobattoIIIC6)
%2017 Maria Lopez-Fernandez

if (RK==1)%BDF1
    s = 1;    
    intflag ='RadauIIA1';    
elseif (RK==2)%Radau IIA of order 3 stage order 2    
    s = 2;    
    intflag ='RadauIIA3';
elseif (RK==3)%Radau IIA order 5, stage order 3    
    s = 3;    
    intflag ='RadauIIA5';
end

Nmax = 3*N;                       % Nmax Power of 2 bigger than 2*N;
Nmax = 2^ceil(log(Nmax)/log(2));

rho = (1e-15)^(1/(2*Nmax));               		 
kv = 0:Nmax-1;                           
xi = rho*exp(1i*2*pi*kv/Nmax);

omegalong = zeros(s,s,Nmax);
%omega = zeros(s,s,N);

for ll=1:Nmax
    xx = xi(ll);
    DM = eval(sprintf('Delta%s(xx)',intflag));
    [V,D] = eig(DM);
    temp = K(diag(D/dt));    
    omegalong(:,:,ll) = V*diag(temp)*inv(V); %K(DM/dt); %LT of the convolution kernel evaluated at \Delta(\zeta/dt).
end

omega0 = zeros(s,s);
omega = zeros(s,s,N);
for kk =1:s
    for jj = 1:s
        omegaV = fft(omegalong(kk,jj,:)); 
        omegaV = squeeze(omegaV);
        omegaV = (rho.^(-kv(1:N)).').*omegaV(1:N);
        omega0(kk,jj) = omegaV(1)/Nmax;
        omega(kk,jj,:)=omegaV/Nmax;
    end
end
DO = eval(sprintf('Delta%s(0)',intflag));
[V,D] = eig(DO);
temp = K(diag(D/dt));
omega0 = V*diag(temp)*inv(V); %K(DM/dt); %LT of the convolution kernel evaluated at \Delta(\zeta/dt).
