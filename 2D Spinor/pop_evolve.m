evo_time = 0.5;
Deltat_count = 1e-4;
Deltat = 5e-5;
evoN = int16(evo_time/Deltat_count);
% phi_save = zeros(Nx,Ny,3,evoN);
fp1 = zeros(evoN,2,4);
for j = 7:8
phi_evo(:,:,1) = phibig(:,:,j);
phi_evo(:,:,2) = zeros(Nx,Ny);
phi_evo(:,:,3) = zeros(Nx,Ny);
%population11 = zeros(3,evoN);
for i = 1:evoN
    phi_evo = dynamic(phi_evo,Deltat_count,Deltat,1*c0,1*c2,Nx,Ny,V,k_scale,fx,deltax,deltay,deltafx,deltafy,Lx,Omega,paritx(1:Nx-1,1:Ny-1,:),parity(1:Nx-1,1:Ny-1,:),dispersion,TF_radius,detuning);
    phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
    
    ff = fourier_transform2(phi_evo,paritx,parity,deltax,deltay);
    ff = ff/norm2d(ff,Nx,Ny,deltafx,deltafy);
    fp1(i,1,j-6) = integr2d(sq(ff(128-20:128+20,128-20:128+20,1)),41,41,deltafx,deltafy);
    fp1(i,2,j-6) = integr2d(sq(ff(128-20:128+20,57-20:57+20,2)),41,41,deltafx,deltafy);
%     phi_save(:,:,:,i) = phi_evo;
    
%     fp00(i,:) = find_peak2d(phi_evo,k_spacing,Nx,Ny,deltafx,deltafy,deltax,deltay,paritx,parity);
    
%     population11(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
%     population11(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
%     population11(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
%    
   
end
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
