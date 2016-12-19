Deltat = 1e-9;
Stop_crit = 1e-9;



potential = @(x)(1/2*x.^2);



phi_0 = Thomas_fermi1D(Beta,X,TF_radius,potential(X),Nx,deltax);
phi = phi_0;
difference = 1;
evo = 200;
n = 0;
while (difference)
    %phi_up = strang_evolve(phi, potential,Deltat,X,Beta,Nx,deltax,deltaf,L );
    phi_up = time_evolve(phi, potential,Deltat,X,Beta,Nx,deltax,deltaf,L );
    if (max(abs(phi_up-phi)) < Stop_crit)
        difference = 0;
    end
    
    if (mod(n,evo) == 0)
        max(abs(phi_up-phi))
        chem_pot(phi,X,Nx,Beta,k_scale,deltax,deltaf,V,L)
    end
    phi = phi_up;
    n = n + 1;
    
end

phi_1 = phi;





