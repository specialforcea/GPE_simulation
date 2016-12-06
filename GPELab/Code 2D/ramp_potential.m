function P = ramp_potential(t,x,y)

t_ramp = 88;

v = 100;
kappa = v/t_ramp;
k = 10;

P = (kappa*t*sin(k*x) .^2);

end

