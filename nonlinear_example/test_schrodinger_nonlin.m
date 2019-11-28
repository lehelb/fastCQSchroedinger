addpath ../fast_schroedinger_codes/
% this example requires error function of complex
% arguments. Here we assume that the Matlab package
% of Marcel Leutenegger is installed.
addpath ../ErrorFunction/
%%Data
sigma=0.8;  
T = 200; %final time
dt = 1/8; %time-step
N = ceil(T/dt); %number of time-steps

K = @(z,d) (exp(1i*pi/4)/2)*z.^(-1/2).*exp(-d*exp(-1i*pi/4)*sqrt(z));
G = @(z,d) 1i*exp(1i*pi/4)*z.^(-1/2).*cosh(d*exp(1i*pi/4)*sqrt(z)); 
x1=-3; x2=3; a=x2; alfa=sqrt(.01); beta=sqrt(.99);
gamma1=-.5; gamma2=-.5; d=2*a;
[lf,le,Nf,Ne]=cnorm(gamma1,gamma2,a);

%For the correction term. Primitive of the convolution kernel
fcorr=@(t,d) exp(1i*pi/4)/(2*sqrt(pi))*(2*sqrt(t).*exp(1i*d^2./(4*t)) + ...
    exp(1i*pi*3/4)*sqrt(pi)*d*erfz(exp(1i*pi*3/4)*d./sqrt(4*t)) - exp(-1i*pi/4)*sqrt(pi)*d);
fcorrd0=@(t) exp(1i*pi/4)*sqrt(t/pi);

RK = 2;
[A,b,c,intflag] = RKdata(RK);
s = length(b);

%The same quadrature can be used for every 0<= d <=2.
B = 3; tol=1e-6;
%Quadrature nodes and weights for the inverse Laplace transform
[Xq,Wq,n0]=quadrature_cqw(d,tol,T,dt,B,RK);
NQ=length(Xq);
%%
%Prepare rs, Ss and values of G
r=eval(sprintf('r_%s(-dt*Xq)',intflag));
S=zeros(NQ,s,s);
for ll=1:NQ
    S(ll,:,:)=eval(sprintf('S_%s(-dt*Xq(ll))',intflag));
end
Gd0=G(Xq,0);
Gd2=G(Xq,d);

q1=zeros(s,N+1); q2=zeros(s,N+1);
tvec=zeros(s,N+1);
[solminusa,solplusa]=freesol(tvec(:,1),a,lf,le,Nf,Ne,alfa,beta);
iniminusa=solminusa(end);
iniplusa=solplusa(end);

q1(:,1)=solminusa;  q2(:,1)=solplusa;
Qhist=zeros(2*s,NQ);
po=2*sigma; 
gamma=2*gamma1/(abs(solminusa(end))^po + abs(solplusa(end))^po);

%For the direct weights
[Wd0,Wd00] = convw_rk(n0+1,dt,RK,@(z) K(z,0));
[Wd2,Wd20] = convw_rk(n0+1,dt,RK,@(z) K(z,d));
Wd0=gamma*Wd0; Wd00=gamma*Wd00;
Wd2=gamma*Wd2; Wd20=gamma*Wd20;

niterv=zeros(1,N+1); 
if po==0
    itermax=1;
else
    itermax=10;
end

h = waitbar(0,'computing');
for n=0:N
    tvec(:,n+1)=n*dt+c*dt;   
       
    [solminusa,solplusa]=freesol(tvec(:,n+1),a,lf,le,Nf,Ne,alfa,beta);
    rhs=[solminusa; solplusa];
    rhs=rhs-gamma*[abs(iniminusa).^po*iniminusa*fcorrd0(tvec(:,n+1)) + abs(iniplusa).^po*iniplusa*fcorr(tvec(:,n+1),d) ;
                   abs(iniplusa).^po*iniplusa*fcorrd0(tvec(:,n+1)) + abs(iniminusa).^po*iniminusa*fcorr(tvec(:,n+1),d)] ;
    rhs=rhs+[Wd00*(abs(iniminusa)^po*iniminusa)*ones(s,1) + Wd20*(abs(iniplusa)^po*iniplusa)*ones(s,1); 
              Wd00*(abs(iniplusa)^po*iniplusa)*ones(s,1) + Wd20*(abs(iniminusa)^po*iniminusa)*ones(s,1)]; 
    hist=0; 
    for ll=1:min(n,n0) %n
        hist = hist + [Wd0(:,:,ll+1)  Wd2(:,:,ll+1);
                      Wd2(:,:,ll+1)  Wd0(:,:,ll+1) ]*...
                [(abs(q1(:,n-ll+1)).^po).*q1(:,n-ll+1)-abs(iniminusa)^po*iniminusa;
                 (abs(q2(:,n-ll+1)).^po).*q2(:,n-ll+1)-abs(iniplusa)^po*iniplusa];
    end
    if n>n0
        for kk=1:NQ
            sqS = squeeze(S(kk,:,:));
            Qhist(:,kk)=r(kk)*Qhist(:,kk)+ [Gd0(kk)*sqS  Gd2(kk)*sqS;
                                            Gd2(kk)*sqS  Gd0(kk)*sqS]*...
                                [(abs(q1(:,n-n0)).^po).*q1(:,n-n0)-abs(iniminusa)^po*iniminusa;
                                 (abs(q2(:,n-n0)).^po).*q2(:,n-n0)-abs(iniplusa)^po*iniplusa];
            hist=hist + gamma*dt*Wq(kk)*r(kk)^(n0+1)*Qhist(:,kk);
        end
    end
    rhs=rhs-hist;
    
    if n==0
        qnew=[iniminusa*ones(s,1); iniplusa*ones(s,1)];
    else
        qnew=[q1(:,n); q2(:,n) ];
    end
    niter=0; err=1;max_err = 0;
    %a simple fixed point iteration
    while ((err>1e-6) && (niter<itermax))
        %System matrix
        M=[eye(s)+Wd00*diag(abs(qnew(1:s)).^po)     Wd20*diag(abs(qnew(s+1:end)).^po)  ;
            Wd20*diag(abs(qnew(1:s)).^po)        eye(s)+Wd00*diag(abs(qnew(s+1:end)).^po)];                
        temp=M\rhs;
        q1(:,n+1)=temp(1:s);
        q2(:,n+1)=temp(s+1:end);
        err=norm(qnew-[q1(:,n+1);q2(:,n+1)],inf);
        qnew = temp;
        niter=niter+1;
        max_err = max(err, max_err);
    end
    niterv(n+1)=niter;        
    waitbar(n/N,h);
end
close(h)
%%
nplot=N+1;figure(1);clf
plot((1:nplot)*dt,abs(q1(s,1:nplot)).^2,'linewidth',2);hold on
plot((1:nplot)*dt,abs(q2(s,1:nplot)).^2,'--','linewidth',2);hold off
a = axis; a(2) = 200; axis(a);
hp = xlabel('$t$');
set(hp,'FontSize',14,'Interpreter','Latex')
hp = ylabel('$|q_j(t)|^2$');set(hp,'FontSize',14,'Interpreter','Latex')
hp = legend('$|q_1(t)|^2$', '$|q_2(t)|^2$'); set(hp,'FontSize',14,'Interpreter','Latex')