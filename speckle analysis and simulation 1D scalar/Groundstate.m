Deltat = 1e-3;
Stop_crit = 3e-16;



potential = @(x)(1/2*x.^2);


init_spin = [1 0 0];
phi_0 = Thomas_fermi1D(c0,X,TF_radius,potential(X),Nx,deltax);
phi = phi_0;
diff = 1;

%show difference every evo of evolutions
evo = 200;
n = 0;
draw = 0;




while (diff)
   
    %phi_up = strang_evolve(phi, potential,Deltat,X,Beta,Nx,deltax,deltaf,L );
    phi_up = time_evolve_ground(phi, potential,Deltat,X,Nx,deltax,deltaf,L,c0);
 
    difference = max(abs(phi_up-phi));
    if (difference < Stop_crit)
        diff = 0;
    end
    
    if (draw && mod(n,evo) == 0)
        plot(X,sq(phi(1,:)))
        drawnow;
%         max(abs(phi_up-phi))
%         chem_pot(phi,X,Nx,c0,k_scale,deltax,deltaf,V,L)
    end
    phi = phi_up;
    n = n + 1;
    
    
end



phi_0 = [init_spin(1).*phi;init_spin(2).*phi;init_spin(3).*phi];




