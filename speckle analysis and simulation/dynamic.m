function phi_2= dynamic(phi,t_evo,Deltat,c0,c2,Nx,speckle,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax)

 
Stop_time = t_evo;










potential = @(x)(0*x.^2);
%order = 2;

X = linspace(xmin,xmax,Nx);
%N_tf = int16(TF_radius/DeltaX);
%cutX = X < TF_radius;
%chem_potential = zeros(1,int16(Stop_time/Deltat));
t = 0;
fftNx = Nx -1;
fftphi = phi(:,1:fftNx);
fftX = X(1:fftNx);
fftL = L - deltax;
fftf = f(1:fftNx);
%speckle = speckle(1:fftNx); 
%fftphi = phi;
Deltat = 1i*Deltat;
evo = 5000;
n = 0;
draw = 1;
width_0 = zeros(1,500);
T = zeros(1,500);
m=1;
while (t < Stop_time)
    %fftphi = strang_evolve(fftphi, potential,Deltat,fftX,Beta,fftNx,deltax,deltaf,fftL);
    fftphi = time_evolve(fftphi, potential,speckle,Deltat,fftX,fftNx,deltax,deltaf,fftL,c0,c2,Omega,k_scale);
    t = t - 1i*Deltat;
    n = n + 1;
    %chem = chem_pot(fftphi,fftX,fftNx,Beta,k_scale,deltax,deltaf,V,fftL);
   %chem_potential(n) = chem;
    
    if (mod(n,evo)==0 && draw == 1)
        width_0(m) = sqrt(integr(fftX.^2.*sq(fftphi(1,:)),fftNx,deltax));
        T(m) = m*-1i*Deltat*evo;
        m = m + 1;
        %fftft = fourier_transform(fftphi,fftNx,deltax);
        %integr(fftf.^2.*sq(fftft(1,:)),fftNx,deltaf)
        %plot(fftf,sq(fftft))
        plot(fftX,sq(fftphi))
        %chem_pot(fftphi(1,:),fftX,fftNx,c0,k_scale,deltax,deltaf,0,fftL)
        
%         fphi = fourier_transform(fftphi,fftNx,deltax);
%         %fphi = fphi./norm1d(fphi,Nx,DeltaX);
%         b = fphi.*conj(fphi);
%         
        %kin = ave_kin(fftphi,fftNx,deltax,deltaf,fftL)
        %kin2 = ave_kin2(fftphi,fftNx,deltax,deltaf,fftL)
         %pot = ave_pot(fftphi,fftNx,deltax,fftX,0,k_scale,c0);
%         plot(fftf,b);
%         title(strcat(num2str(chem),'||',num2str(kin),'||',num2str(pot),'||',num2str(kin2)));
        drawnow;
        
    end
    
end
phi_2 = zeros(3,Nx);
phi_2(:,1:fftNx) = fftphi;
phi_2(:,Nx) = phi_2(:,1);
save('width_0.mat','width_0');
save('T0.mat','T');

end