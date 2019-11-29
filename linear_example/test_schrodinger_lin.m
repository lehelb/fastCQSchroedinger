%%
addpath ../fast_schroedinger_codes/

%Data
V1bar=-(1+tanh(1)); V2bar=V1bar; 
V1=@(t) V1bar*((1+sin(t)).*heaviside(t)+heaviside(-t-1e-16));
V2=@(t) V2bar*((1-sin(t)).*heaviside(t)+heaviside(-t-1e-16));
psi0 = @(x) (cosh(x)/cosh(1)).*heaviside(1-abs(x)-1e-16)+...
            exp(1-abs(x)).*heaviside(abs(x)-1);

K = @(z,d) .5*exp(1i*pi/4)*z.^(-1/2).*exp(-d*exp(-1i*pi/4)*sqrt(z));
G = @(z,d) 1i*exp(1i*pi/4)*z.^(-1/2).*cosh(d*exp(1i*pi/4)*sqrt(z)); 
x1=-1; x2=1; omega=1;

RK = 2;
[A,b,c] = RKdata(RK);
s = length(b);

T = 40;
dt = .1; N = ceil(T/dt);
d = 2; %The same quadrature can be used for every 0<= d <=2.
B = 3;
tol=1e-6;
%For the new quadrature
[Xq,Wq,n0,Q1,Q2s]=quadrature_cqw(d,tol,T,dt,B,RK);
NQ=length(Xq);
%Prepare rs, Ss and values of G
r = rRadauIIA(-dt*Xq,RK);
S=zeros(NQ,s,s);
for ll=1:NQ
    S(ll,:,:) = SRadauIIA(-dt*Xq(ll),RK);
end
Gd0=G(Xq,0);
Gd2=G(Xq,2);
%For the direct weights
[Wd0,Wd00] = convw_rk(n0+1,dt,RK,@(z) K(z,0));
[Wd2,Wd20] = convw_rk(n0+1,dt,RK,@(z) K(z,d));

q1=zeros(s,N+1); q2=zeros(s,N+1);
tvec=zeros(s,N+1);
tvec(:,1)=-dt+c*dt;
q1(:,1)=psi0(x1)*exp(1i*omega*tvec(:,1));
q2(:,1)=psi0(x2)*exp(1i*omega*tvec(:,1));
Qhist=zeros(2*s,NQ);

wb = waitbar(0,'computing');
for n=1:N
    tvec(:,n+1)=tvec(end,n)+c*dt;
    %System matrix
    M=[eye(s)+Wd00*diag(V1(tvec(:,n+1)))     Wd20*diag(V2(tvec(:,n+1)))  ;
        Wd20*diag(V1(tvec(:,n+1)))        eye(s)+Wd00*diag(V2(tvec(:,n+1)))];
      
    rhs=[(eye(s)+Wd0(:,:,1)*V1bar)*psi0(x1)*exp(1i*omega*tvec(:,n+1)) + Wd2(:,:,1)*V2bar*psi0(x2)*exp(1i*omega*tvec(:,n+1))  ;
         (eye(s)+Wd0(:,:,1)*V2bar)*psi0(x2)*exp(1i*omega*tvec(:,n+1)) + Wd2(:,:,1)*V1bar*psi0(x1)*exp(1i*omega*tvec(:,n+1))];
    
    hist=0;     
    for ll=1:min(n,n0) 
        hist = hist + [Wd0(:,:,ll+1)  Wd2(:,:,ll+1);
            Wd2(:,:,ll+1)  Wd0(:,:,ll+1) ]*...
            [V1(tvec(:,n-ll+1)).*q1(:,n-ll+1)-V1bar*psi0(x1)*exp(1i*omega*(tvec(:,n-ll+1)));
             V2(tvec(:,n-ll+1)).*q2(:,n-ll+1)-V2bar*psi0(x2)*exp(1i*omega*(tvec(:,n-ll+1)))];
    end
    if n>n0
        for kk=1:NQ
            sqS = squeeze(S(kk,:,:));
            Qhist(:,kk)=r(kk)*Qhist(:,kk)+ [Gd0(kk)*sqS  Gd2(kk)*sqS;
                                            Gd2(kk)*sqS  Gd0(kk)*sqS]*...
                         [V1(tvec(:,n-n0)).*q1(:,n-n0) - V1bar*psi0(x1)*exp(1i*omega*(tvec(:,n-n0)));
                          V2(tvec(:,n-n0)).*q2(:,n-n0) - V2bar*psi0(x2)*exp(1i*omega*(tvec(:,n-n0)))];
            hist=hist + dt*Wq(kk)*r(kk)^(n0+1)*Qhist(:,kk);
        end
    end
    
    rhs=rhs-hist;
    temp=M\rhs;    
    q1(:,n+1)=temp(1:s);
    q2(:,n+1)=temp(s+1:end);
    waitbar(n/N,wb);
end
close(wb)
nplot=N+1;

figure(1);clf
plot((1:nplot)*dt,abs(q1(s,1:nplot)),'linewidth',2); hold on
plot((1:nplot)*dt,abs(q2(s,1:nplot)),'--','linewidth',2); hold off
hp = xlabel('$t$');
set(hp,'FontSize',14,'Interpreter','Latex')
hp = ylabel('$|q_j(t)|$');set(hp,'FontSize',14,'Interpreter','Latex')
a = axis; a(2) = T; axis(a);
hp = legend('$|q_1(t)|$', '$|q_2(t)|$'); set(hp,'FontSize',14,'Interpreter','Latex')