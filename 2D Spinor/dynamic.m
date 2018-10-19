function phi_2 = dynamic(phi,t_evo,Deltat,c0,c2,X,Y,Nx,Ny,trap,V,k_scale,deltax,deltay,deltafx,deltafy,fx,fy,Omega,paritx,parity,dispersion,detuning)

 
Stop_time = t_evo;










potential = @(x,y)(trap/2*x.^2 + trap/2*y.^2);
%order = 2;



t = 0;
fftNx = Nx -1;
fftNy = Ny -1;
fftphi = phi(1:fftNx,1:fftNy,:);
fftX = X(1:fftNx,1:fftNy);
fftY = Y(1:fftNx,1:fftNy);

fftfx = fx(1:fftNx);
fftfy = fy(1:fftNx);
Deltat = 1i*Deltat;
evo = 1;
n = 0;
draw = 0;

while (t < Stop_time)
    
    fftphi = time_evolve(fftphi, potential,Deltat,fftX,fftY,fftNx,fftNy,deltax,deltay,deltafx,deltafy,c0,c2,Omega,k_scale,paritx,parity,dispersion,detuning);
    t = t - 1i*Deltat;
    n = n + 1;
  
    
    if (mod(n,evo)==0 && draw == 1)
        fphi = fourier_transform2(fftphi,paritx,parity,deltax,deltay);
       
        b = fphi.*conj(fphi);
        

       
        
        
    end
    
end
phi_2 = zeros(Nx,Ny,3);
phi_2(1:fftNx,1:fftNy,:) = fftphi;
phi_2(Nx,:,:) = phi_2(1,:,:);
phi_2(:,Ny,:) = phi_2(:,1,:);

end