function gamma = test_gamma(xis,RK)

[A,b,~] = RKdata(RK);
m = length(b);    
gamma = zeros(size(xis));
r = @(z) 1+z*b'*((eye(m)-z*A)\ones(m,1));

for k = 1:length(xis)
    xi = xis(k);
    z = -xi+linspace(0,10,10)*1i;
    res = zeros(length(z),1);
    for j = 1:length(z)
        res(j) = abs(r(z(j)));
    end
    gamma(k) = log(1/max(res))/xi;
end