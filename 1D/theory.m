t = (1:1:200).*10^-6;
alpha = E_r/hbar.*t;
V_0 = 4*E_r;
beta = V_0/hbar.*t;
%plot(beta,besselj(0,beta./2).^2,beta,besselj(1,beta./2).^2,beta,besselj(2,beta./2).^2,beta,besselj(3,beta./2).^2,beta,besselj(4,beta./2).^2)
plot(beta,besselj(20,beta./2).^2)