function S = S_RadauIIA5(z) 

% S(z) = r^-1 * (I - z*A) *e* bT * (I - z*A)   

S = ([[-2./9.*(9+6.^(1./2)).*(3.*z+5.*6.^(1./2)).*(-10-4.*z+z.^2-z.*6.^(1./2)+10.*6.^(1./2))./(20+8.*z+z.^2)./(-60+36.*z-9.*z.^2+z.^3),2./45.*(51+11.*6.^(1./2)).*(3.*z-5.*6.^(1./2)).*(-10-4.*z+z.^2-z.*6.^(1./2)+10.*6.^(1./2))./(20+8.*z+z.^2)./(-60+36.*z-9.*z.^2+z.^3),-2./15.*(1+6.^(1./2)).*(3.*z.^2-12.*z+20).*(-10-4.*z+z.^2-z.*6.^(1./2)+10.*6.^(1./2))./(20+8.*z+z.^2)./(-60+36.*z-9.*z.^2+z.^3)];
    [-2./45.*(-51+11.*6.^(1./2)).*(-4.*z+z.^2+z.*6.^(1./2)-10-10.*6.^(1./2)).*(3.*z+5.*6.^(1./2))./(20+8.*z+z.^2)./(-60+36.*z-9.*z.^2+z.^3),2./9.*(-9+6.^(1./2)).*(-4.*z+z.^2+z.*6.^(1./2)-10-10.*6.^(1./2)).*(3.*z-5.*6.^(1./2))./(20+8.*z+z.^2)./(-60+36.*z-9.*z.^2+z.^3),2./15.*(-1+6.^(1./2)).*(-4.*z+z.^2+z.*6.^(1./2)-10-10.*6.^(1./2)).*(3.*z.^2-12.*z+20)./(20+8.*z+z.^2)./(-60+36.*z-9.*z.^2+z.^3)];
    [-1./9.*(-3+8.*6.^(1./2)).*(3.*z+5.*6.^(1./2))./(-60+36.*z-9.*z.^2+z.^3),1./9.*(3+8.*6.^(1./2)).*(3.*z-5.*6.^(1./2))./(-60+36.*z-9.*z.^2+z.^3),-1./3.*(3.*z.^2-12.*z+20)./(-60+36.*z-9.*z.^2+z.^3)]]);
