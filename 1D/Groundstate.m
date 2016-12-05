Deltat = 1e-4;
Stop_crit = 1e-8;

Beta = 3000; 
TF_radius = (3*Beta/2)^(1/3);
xmin = -TF_radius;
xmax = TF_radius;

Nx = 2^11+1;


Delta = 0.5;


DeltaX = (xmax-xmin)/(Nx-1);

potential = @(x)(1/2*x.^2);


X = linspace(xmin,xmax,Nx);
N_tf = int16(TF_radius/DeltaX);
cutX = X < TF_radius;
phi_0 = Thomas_fermi1D(Beta,potential(X)).*cutX;

phi = phi_0;
difference = 1;
while (difference)
    phi_up = time_evolve(phi, potential,Deltat,X,Beta,Nx,DeltaX);
    if (max(abs(phi_up-phi)) < Stop_crit)
        difference = 0;
    end
    max(abs(phi_up-phi))
    phi = phi_up;
    
end

phi_1 = phi;






