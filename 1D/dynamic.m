Deltat = 1e-5;
Stop_time = 0.05;

Beta = 3000; 
TF_radius = (3*Beta/2)^(1/3);
xmin = -TF_radius;
xmax = TF_radius;

Nx = 2^11+1;


Delta = 0.5;

lattice_pot = 4*E_r;
V = lattice_pot/quanta;
DeltaX = (xmax-xmin)/(Nx-1);
k_scale = k*a_0;
potential = @(x)(1/2*x.^2 + V*sin(k_scale*x).^2 );


X = linspace(xmin,xmax,Nx);
N_tf = int16(TF_radius/DeltaX);
cutX = X < TF_radius;

t = 0;
phi = phi_1;
Deltat = 1i*Deltat;
while (t < Stop_time)
    phi = time_evolve(phi, potential,Deltat,X,Beta,Nx,DeltaX);
    t = t - 1i*Deltat;
    t
    
    
end
phi_2 = phi;