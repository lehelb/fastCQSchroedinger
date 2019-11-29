function [omega,omega0] = convw_rk(N,dt,RK,K)
%function f = convw_rk(N,dt,RK,Ks)
%returns conv. weight \omega_j(K) for j = 0, .., N
%RK = (1 RadauIIA1) (2 RadauIIA2) (3 RadauIIA3) (4 LobattoIIIC6)
%2017 Maria Lopez-Fernandez

s = RK; %number of stages

Nmax = 3*N;                       % Nmax Power of 2 bigger than 2*N;
Nmax = 2^ceil(log(Nmax)/log(2));

rho = (1e-15)^(1/(2*Nmax));               		 
kv = 0:Nmax-1;                           
xi = rho*exp(1i*2*pi*kv/Nmax);

omegalong = zeros(s,s,Nmax);
%omega = zeros(s,s,N);
for ll=1:Nmax
    xx = xi(ll);
    DM = deltaRadauIIA(xx,RK);
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
DO = deltaRadauIIA(0,RK);
[V,D] = eig(DO);
temp = K(diag(D/dt));
omega0 = V*diag(temp)*inv(V); %K(DM/dt); %LT of the convolution kernel evaluated at \Delta(\zeta/dt).