Deltat = 1e-9;
Stop_crit = 2.5e-8;



potential = @(x)(1/2*x.^2);



phi_0 = Thomas_fermi1D(c0,X,TF_radius,potential(X),Nx,deltax);
phi = phi_0;
diff = 1;

%show some results every evo of evolutions, difference is shown after every
%single evolution
evo = 200;
n = 0;





while (diff)
   
    
    phi_up = time_evolve_ground(phi, potential,Deltat,X,Nx,deltax,deltaf,L,c0);
    difference = max(abs(phi_up-phi))
    if (difference < Stop_crit)
        diff = 0;
    end
    
%     if (mod(n,evo) == 0)
%         max(abs(phi_up-phi))
%         chem_pot(phi,X,Nx,c0,k_scale,deltax,deltaf,V,L)
%     end
    phi = phi_up;
    n = n + 1;
    
    
end



phi_0 = phi;





