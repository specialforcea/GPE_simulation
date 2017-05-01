evo_time = 0.3;
% time = [0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1];
Deltat_count = 1e-4;
% timeN = time./Deltat_count;
Deltat = 1e-5;
evoN = int16(evo_time/Deltat_count);
%phifindout = zeros(Nx,Ny,3,evoN);
%phi_half = zeros(Nx,Ny,3,evoN);
% ksp = zeros(Nx,Ny,2,50,4);

%fp1 = zeros(evoN,2,4);
real_Omega = [1 2 3 4 5 6]*2.34*1e4;
det = [1 2 4 5 50];
%fpr = zeros(evoN,3,6,5);
%for j = 1:3

% nn = 0;
Omega = real_Omega/Dip_freq;
for yy = 1:6
for yyy = 2:5
    phi_evo = phi;

for i = 1:evoN
    phi_evo = dynamic(phi_evo,Deltat_count,Deltat,1*c0,1*c2,Nx,Ny,V,k_scale,fx,deltax,deltay,deltafx,deltafy,Lx,Omega(yy),paritx(1:Nx-1,1:Ny-1,:),parity(1:Nx-1,1:Ny-1,:),dispersion,TF_radius,det(yyy)*detuning);
    phi_evo = phi_evo./norm2d(phi_evo, Nx,Ny, deltax,deltay);
    
    
    
    %phifindout(:,:,:,i) = phi_evo;
    

    
%     if mod(i,100) == 0
%         nn = nn + 1;
%         ff = fourier_transform2(phi_evo,paritx,parity,deltax,deltay);
%         ff = ff/norm2d(ff,Nx,Ny,deltafx,deltafy);
%         fpr(i,1,j-6) = integr2d(sq(ff(128-20:128+20,128-20:128+20,1)),41,41,deltafx,deltafy);
%         fpr(i,2,j-6) = integr2d(sq(ff(128-20:128+20,57-20:57+20,2)),41,41,deltafx,deltafy);
% %         ksp(:,:,1,nn,j) = ff(:,:,1);
% %         ksp(:,:,2,nn,j) = ff(:,:,2);
%     end
    ff = fourier_transform2(phi_evo,paritx,parity,deltax,deltay);
    ff = ff/norm2d(ff,Nx,Ny,deltafx,deltafy);
    [colmax,rowind] = max(sq(ff(:,:,1)));
    [~,colind] = max(colmax);
    fpr(i,1,yy,yyy) = integr2d(sq(ff(max(1,rowind(colind)-20):min(Ny,rowind(colind)+20),max(1,colind-20):min(Nx,colind+20),1)),size(max(1,rowind(colind)-20):min(Ny,rowind(colind)+20),2),size(max(1,colind-20):min(Nx,colind+20),2),deltafx,deltafy);
    [colmax,rowind] = max(sq(ff(:,:,2)));
    [~,colind] = max(colmax);
    fpr(i,2,yy,yyy) = integr2d(sq(ff(max(1,rowind(colind)-20):min(Ny,rowind(colind)+20),max(1,colind-20):min(Nx,colind+20),2)),size(max(1,rowind(colind)-20):min(Ny,rowind(colind)+20),2),size(max(1,colind-20):min(Nx,colind+20),2),deltafx,deltafy);
    [colmax,rowind] = max(sq(ff(:,:,3)));
    [~,colind] = max(colmax);
    fpr(i,3,yy,yyy) = integr2d(sq(ff(max(1,rowind(colind)-20):min(Ny,rowind(colind)+20),max(1,colind-20):min(Nx,colind+20),3)),size(max(1,rowind(colind)-20):min(Ny,rowind(colind)+20),2),size(max(1,colind-20):min(Nx,colind+20),2),deltafx,deltafy);
%     phi_save(:,:,:,i) = phi_evo;
    
%     fp00(i,:) = find_peak2d(phi_evo,k_spacing,Nx,Ny,deltafx,deltafy,deltax,deltay,paritx,parity);
    
%     population11(1,i) = integr2d(sq(phi_evo(:,:,1)),Nx,Ny,deltax,deltay);
%     population11(2,i) = integr2d(sq(phi_evo(:,:,2)),Nx,Ny,deltax,deltay);
%     population11(3,i) = integr2d(sq(phi_evo(:,:,3)),Nx,Ny,deltax,deltay);
%    
   
end
end
end
save('fpr.mat','fpr');
% phi_half(:,:,:,j-6) = phi_evo;
%end


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
