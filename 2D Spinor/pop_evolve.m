evo_time = 0.12;
Deltat_count = 1e-4;
Deltat = 5e-5;
evoN = int16(evo_time/Deltat_count);
population11 = zeros(3,evoN);
phi_evo = phi_0;
for i = 1:evoN
    phi_evo = dynamic(phi_evo,Deltat_count,Deltat,c0,c2,Nx,Ny,V,k_scale,f,deltax,deltay,deltafx,deltafy,L,Omega,paritx,parity,dispersion);
    phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
    population11(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
    population11(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
    population11(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
    %population(1,i) = sq(phi_evo(100,100,1));
   
end


% popN = linspace(1e-4,1e-4*399,400);
% plot(popN,population(:,1:400)')
population1001 = zeros(3,evoN);
phi_evo = phi_0;
for i = 1:evoN
    phi_evo = dynamic(phi_evo,Deltat_count,Deltat,100*c0,c2,Nx,Ny,V,k_scale,f,deltax,deltay,deltafx,deltafy,L,Omega,paritx,parity,dispersion);
    phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
    population1001(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
    population1001(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
    population1001(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
    %population(1,i) = sq(phi_evo(100,100,1));
   
end

population10001 = zeros(3,evoN);
phi_evo = phi_0;
for i = 1:evoN
    phi_evo = dynamic(phi_evo,Deltat_count,Deltat,1000*c0,c2,Nx,Ny,V,k_scale,f,deltax,deltay,deltafx,deltafy,L,Omega,paritx,parity,dispersion);
    phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
    population10001(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
    population10001(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
    population10001(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
    %population(1,i) = sq(phi_evo(100,100,1));
   
end

population1100 = zeros(3,evoN);
phi_evo = phi_0;
for i = 1:evoN
    phi_evo = dynamic(phi_evo,Deltat_count,Deltat,c0,100*c2,Nx,Ny,V,k_scale,f,deltax,deltay,deltafx,deltafy,L,Omega,paritx,parity,dispersion);
    phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
    population1100(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
    population1100(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
    population1100(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
    %population(1,i) = sq(phi_evo(100,100,1));
   
end

population11000 = zeros(3,evoN);
phi_evo = phi_0;
for i = 1:evoN
    phi_evo = dynamic(phi_evo,Deltat_count,Deltat,c0,1000*c2,Nx,Ny,V,k_scale,f,deltax,deltay,deltafx,deltafy,L,Omega,paritx,parity,dispersion);
    phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
    population11000(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
    population11000(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
    population11000(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
    %population(1,i) = sq(phi_evo(100,100,1));
   
end

population15000 = zeros(3,evoN);
phi_evo = phi_0;
for i = 1:evoN
    phi_evo = dynamic(phi_evo,Deltat_count,Deltat,c0,5000*c2,Nx,Ny,V,k_scale,f,deltax,deltay,deltafx,deltafy,L,Omega,paritx,parity,dispersion);
    phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
    population15000(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
    population15000(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
    population15000(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
    %population(1,i) = sq(phi_evo(100,100,1));
   
end
