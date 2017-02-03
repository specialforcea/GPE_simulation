Deltat = 1e-8;
Stop_crit = 1e-8;



potential = @(x,y)(1/2*x.^2 + 1/2*y.^2);


init_spin = [1 0 0];
phi_0 = Thomas_fermi1D(c0,X,TF_radius,potential(X),Nx,deltax);
phi = phi_0;
difference = 1;
evo = 200;
n = 0;
while (difference)
    %phi_up = strang_evolve(phi, potential,Deltat,X,Beta,Nx,deltax,deltaf,L );
    phi_up = time_evolve(phi, potential,Deltat,X,c0,Nx,deltax,deltaf,L );
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

phi_0 = [init_spin(1).*phi;init_spin(2).*phi;init_spin(3).*phi];




