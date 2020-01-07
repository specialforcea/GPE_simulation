Deltat = 1i*1e-5;
Stop_crit = 3e-16;



potential = @(x)(1/2*x.^2);



phi = phi_0;
diff = 1;

%show difference every evo of evolutions
evo = 200;
n = 0;





while (diff)
   
    %phi_up = strang_evolve(phi, potential,Deltat,X,Beta,Nx,deltax,deltaf,L );
    phi_up = time_evolve(phi, potential, 0,Deltat, X, Nx, deltax,deltaf,L,c0,c2,detuning,k_R,detuning);

    difference = max(max(abs(phi_up-phi)))
    if (difference < Stop_crit)
        diff = 0;
    end
    
    if(mod(n,evo) == 0)
        plot(X,sq(phi))
        drawnow;
    end
    
%     if (mod(n,evo) == 0)
%         max(abs(phi_up-phi))
%         chem_pot(phi,X,Nx,c0,k_scale,deltax,deltaf,V,L)
%     end
    phi = phi_up;
    n = n + 1;
    
    
end
