evo_time = 0.20;
Deltat_count = 0.05;
Deltat = 5e-5;
evoN = int16(evo_time/Deltat_count);
phi_save = zeros(Nx,Ny,3,evoN);
% fp00 = zeros(evoN,2);
phi_evo = phi_0;
for i = 1:evoN
    phi_evo = dynamic(phi_evo,Deltat_count,Deltat,10*c0,1*c2,Nx,Ny,V,k_scale,fx,deltax,deltay,deltafx,deltafy,Lx,Omega,paritx,parity,dispersion,TF_radius,detuning);
    phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
    
    phi_save(:,:,:,i) = phi_evo;
    
%     fp00(i,:) = find_peak2d(phi_evo,k_spacing,Nx,Ny,deltafx,deltafy,deltax,deltay,paritx,parity);
    
%     population00(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
%     population00(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
%     population00(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
    %population(1,i) = sq(phi_evo(100,100,1));
   
end


% popN = linspace(1e-4,1e-4*399,400);
% plot(popN,population(:,1:400)')

% population11 = zeros(3,evoN);
% phi_evo = phi_0;
% for i = 1:evoN
%     phi_evo = dynamic(phi_evo,Deltat_count,Deltat,c0,c2,Nx,Ny,V,k_scale,fx,deltax,deltay,deltafx,deltafy,Lx,Omega,paritx,parity,dispersion,TF_radius,detuning);
%     phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
%     population11(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
%     population11(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
%     population11(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
%     population(1,i) = sq(phi_evo(100,100,1));
%    
% end
% 
% population101 = zeros(3,evoN);
% phi_evo = phi_0;
% for i = 1:evoN
%     phi_evo = dynamic(phi_evo,Deltat_count,Deltat,10*c0,c2,Nx,Ny,V,k_scale,fx,deltax,deltay,deltafx,deltafy,Lx,Omega,paritx,parity,dispersion,TF_radius,detuning);
%     phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
%     population101(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
%     population101(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
%     population101(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
%     population(1,i) = sq(phi_evo(100,100,1));
%    
% end
% 
% population1001 = zeros(3,evoN);
% phi_evo = phi_0;
% for i = 1:evoN
%     phi_evo = dynamic(phi_evo,Deltat_count,Deltat,100*c0,c2,Nx,Ny,V,k_scale,fx,deltax,deltay,deltafx,deltafy,Lx,Omega,paritx,parity,dispersion,TF_radius,detuning);
%     phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
%     population1001(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
%     population1001(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
%     population1001(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
%     population(1,i) = sq(phi_evo(100,100,1));
%    
% end
% 
% 
% 
% population1100 = zeros(3,evoN);
% phi_evo = phi_0;
% for i = 1:evoN
%     phi_evo = dynamic(phi_evo,Deltat_count,Deltat,c0,100*c2,Nx,Ny,V,k_scale,fx,deltax,deltay,deltafx,deltafy,Lx,Omega,paritx,parity,dispersion,TF_radius,detuning);
%     phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
%     population1100(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
%     population1100(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
%     population1100(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
%     population(1,i) = sq(phi_evo(100,100,1));
%    
% end
% 
% population11000 = zeros(3,evoN);
% phi_evo = phi_0;
% for i = 1:evoN
%     phi_evo = dynamic(phi_evo,Deltat_count,Deltat,c0,1000*c2,Nx,Ny,V,k_scale,fx,deltax,deltay,deltafx,deltafy,Lx,Omega,paritx,parity,dispersion,TF_radius,detuning);
%     phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
%     population11000(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
%     population11000(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
%     population11000(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
%     population(1,i) = sq(phi_evo(100,100,1));
%    
% end
% 
% 
