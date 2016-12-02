Deltat = 1e-2;
Stop_crit = 1e-7;

xmin = -2;
xmax = 2;

Nx = 2^10+1;


Delta = 0.5;
Beta = 3000; 

DeltaX = (xmax-xmin)/(Nx-1);

potential = @(x)(1/2*x.^2);

TF_radius = (3*Beta/2)^(1/3);
X = linspace(xmin,xmax,Nx);
N_tf = int16(TF_radius/DeltaX);
cutX = X<TF_radius;
phi_0 = Thomas_fermi1D(Beta,potential(X)).*cutX;





