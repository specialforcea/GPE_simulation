function phi_2= delete_dynamic(phi,t_evo,Deltat,c0,c2,Nx,speckle,trap,V,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax)

 
Stop_time = t_evo;










potential = @(x)(trap/2*x.^2 + V.*sin(k_scale.*x).^2 );
%order = 2;
X = linspace(xmin,xmax,Nx);
%N_tf = int16(TF_radius/DeltaX);
%cutX = X < TF_radius;
%chem_potential = zeros(1,int16(Stop_time/Deltat));
t = 0;
fftNx = Nx;
fftphi = phi;
fftX = X;
fftL = L;
fftf = f;
if speckle~=0
    speckle = speckle(1:fftNx); 
end
%fftphi = phi;
Deltat = 1i*Deltat;
evo = 5000;
n = 0;
draw = 0;
% width_0 = zeros(1,500);
% T = zeros(1,500);
m=1;
while (t < Stop_time)
    %fftphi = strang_evolve(fftphi, potential,Deltat,fftX,Beta,fftNx,deltax,deltaf,fftL);
    fftphi = time_evolve(fftphi, potential,speckle,Deltat,fftX,fftNx,deltax,deltaf,fftL,c0,c2,Omega,k_scale);
    t = t - 1i*Deltat;
    n = n + 1;
    %chem = chem_pot(fftphi(1,:),fftX,fftNx,c0,k_scale,deltax,deltaf,0,fftL)
   %chem_potential(n) = chem;
    
    if (mod(n,evo)==0 && draw == 1)
%         width_0(m) = sqrt(integr(fftX.^2.*sq(fftphi(1,:)),fftNx,deltax));
%         T(m) = m*-1i*Deltat*evo;
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
phi_2 = fftphi;
% save('width_0.mat','width_0');
% save('T0.mat','T');

end