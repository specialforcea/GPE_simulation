function phi_2 = dynamic(phi,t_evo,Deltat,c0,c2,X,Y,Nx,Ny,V,k_scale,f,deltax,deltay,deltafx,deltafy,L,Omega,paritx,parity,dispersion,detuning)

 
Stop_time = t_evo;


%TF_radius = (4*c0/pi)^(1/4);







potential = @(x,y)(1/2*x.^2 + 1/2*y.^2);
%order = 2;


%N_tf = int16(TF_radius/DeltaX);
%cutX = X < TF_radius;
%chem_potential = zeros(1,int16(Stop_time/Deltat));
t = 0;
fftNx = Nx -1;
fftNy = Ny -1;
fftphi = phi(1:fftNx,1:fftNy,:);
fftX = X(1:fftNx,1:fftNy);
fftY = Y(1:fftNx,1:fftNy);
fftL = L - deltax;
fftf = f(1:fftNx);
%fftphi = phi;
Deltat = 1i*Deltat;
evo = 1;
n = 0;
draw = 0;

while (t < Stop_time)
    %fftphi = strang_evolve(fftphi, potential,Deltat,fftX,Beta,fftNx,deltax,deltaf,fftL);
    fftphi = time_evolve(fftphi, potential,Deltat,fftX,fftY,fftNx,fftNy,deltax,deltay,deltafx,deltafy,fftL,c0,c2,Omega,k_scale,paritx,parity,dispersion,detuning);
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
phi_2 = zeros(Nx,Ny,3);
phi_2(1:fftNx,1:fftNy,:) = fftphi;
phi_2(Nx,:,:) = phi_2(1,:,:);
phi_2(:,Ny,:) = phi_2(:,1,:);
%phi_2 = fftphi;
end