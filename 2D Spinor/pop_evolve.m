evo_time = 0.12;
Deltat_count = 1e-4;
Deltat = 5e-5;
evoN = int16(evo_time/Deltat_count);
population00 = zeros(3,evoN);
phi_evo = phi_0;
for i = 1:evoN
    phi_evo = dynamic(phi_evo,Deltat_count,Deltat,0*c0,0*c2,Nx,Ny,V,k_scale,fx,deltax,deltay,deltafx,deltafy,Lx,Omega,paritx,parity,dispersion,TF_radius);
    phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
    population00(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
    population00(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
    population00(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
    %population(1,i) = sq(phi_evo(100,100,1));
   
end


% popN = linspace(1e-4,1e-4*399,400);
% plot(popN,population(:,1:400)')

