function phi_2 = dynamic(phi,t_evo,Deltat,c0,c2,X,Y,Z,Nx,Ny,Nz,V,k_scale,fx,fy,fz,deltax,deltay,deltaz,deltafx,deltafy,deltafz,L,Omega,paritx,parity,paritz,dispersion,TF_radius,detuning)

 
Stop_time = t_evo;





potential = @(x,y,z)(1/2*x.^2 + 1/2*y.^2 + 1/2*z.^2);
%order = 2;

%N_tf = int16(TF_radius/DeltaX);
%cutX = X < TF_radius;
%chem_potential = zeros(1,int16(Stop_time/Deltat));
t = 0;
fftNx = Nx -1;
fftNy = Ny -1;
fftNz = Nz -1;
fftphi = phi(1:fftNy,1:fftNx,1:fftNz,:);
fftX = X(1:fftNy,1:fftNx,1:fftNz);
fftY = Y(1:fftNx,1:fftNy,1:fftNz);
fftZ = Z(1:fftNx,1:fftNy,1:fftNz);

fftL = L - deltax;
fftfx = fx(1:fftNx);
fftfy = fy(1:fftNy);
fftfz = fz(1:fftNz);
%fftphi = phi;
Deltat = 1i*Deltat;
evo = 1;
n = 0;
draw = 0;

while (t < Stop_time)
    %fftphi = strang_evolve(fftphi, potential,Deltat,fftX,Beta,fftNx,deltax,deltaf,fftL);
    fftphi = time_evolve(fftphi, potential,Deltat,fftX,fftY,fftZ,fftNx,fftNy,fftNz,deltax,deltay,deltaz,deltafx,deltafy,deltafz,fftL,c0,c2,Omega,k_scale,paritx,parity,paritz,dispersion,detuning);
    t = t - 1i*Deltat;
    n = n + 1;
    %chem = chem_pot(fftphi,fftX,fftNx,Beta,k_scale,deltax,deltaf,V,fftL);
   %chem_potential(n) = chem;
    
    if (mod(n,evo)==0 && draw == 1)
        fphi = fourier_transform2(fftphi,paritx,parity,deltax,deltay);
        %fphi = fphi./norm1d(fphi,Nx,DeltaX);
        b = fphi.*conj(fphi);
        
%         kin = ave_kin(fftphi,fftNx,deltax,deltaf,fftL);
%         kin2 = ave_kin2(fftphi,fftNx,deltax,deltaf,fftL);
%         pot = ave_pot(fftphi,fftNx,deltax,fftX,V,k_scale,Beta);
       % surf(fftf,fftf,b(:,:,2));
       % title(strcat(num2str(chem),'||',num2str(kin),'||',num2str(pot),'||',num2str(kin2)));
       plot(fftf,b(fftNx/2,:,2)); 
       drawnow;
        
        
    end
    
end
phi_2 = zeros(Ny,Nx,Nz,3);
phi_2(1:fftNy,1:fftNx,1:fftNz,:) = fftphi;
phi_2(Ny,:,:,:) = phi_2(1,:,:,:);
phi_2(:,Nx,:,:) = phi_2(:,1,:,:);
phi_2(:,:,Nz,:) = phi_2(:,:,1,:);
%phi_2 = fftphi;
end