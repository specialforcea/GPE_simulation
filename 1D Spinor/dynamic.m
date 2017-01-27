function phi_2 = dynamic(phi,t_evo,Deltat,c0,c2,Nx,V,k_scale,f,deltax,deltaf,L)

 
Stop_time = t_evo;


TF_radius = (3*c0/2)^(1/3);
xmin = -TF_radius;
xmax = TF_radius;







potential = @(x)(1/2*x.^2);
%order = 2;

X = linspace(xmin,xmax,Nx);
%N_tf = int16(TF_radius/DeltaX);
%cutX = X < TF_radius;
%chem_potential = zeros(1,int16(Stop_time/Deltat));
t = 0;
fftNx = Nx -1;
fftphi = phi(1:fftNx);
fftX = X(1:fftNx);
 fftL = L - deltax;
 fftf = f(1:fftNx);
%fftphi = phi;
Deltat = 1i*Deltat;
evo = 500;
n = 0;
draw = 0;

while (t < Stop_time)
    %fftphi = strang_evolve(fftphi, potential,Deltat,fftX,Beta,fftNx,deltax,deltaf,fftL);
    fftphi = time_evolve(fftphi, potential,Deltat,fftX,Beta,fftNx,deltax,deltaf,fftL);
    t = t - 1i*Deltat;
    n = n + 1;
    chem = chem_pot(fftphi,fftX,fftNx,Beta,k_scale,deltax,deltaf,V,fftL);
   %chem_potential(n) = chem;
    
    if (mod(n,evo)==0 && draw == 1)
        fphi = fourier_transform(fftphi,fftNx,deltax);
        %fphi = fphi./norm1d(fphi,Nx,DeltaX);
        b = fphi.*conj(fphi);
        
        kin = ave_kin(fftphi,fftNx,deltax,deltaf,fftL);
        kin2 = ave_kin2(fftphi,fftNx,deltax,deltaf,fftL);
        pot = ave_pot(fftphi,fftNx,deltax,fftX,V,k_scale,Beta);
        plot(fftf,b);
        title(strcat(num2str(chem),'||',num2str(kin),'||',num2str(pot),'||',num2str(kin2)));
        drawnow;
        
    end
    
end
phi_2 = zeros(1,Nx);
phi_2(1:fftNx) = fftphi;
phi_2(Nx) = phi_2(1);
%phi_2 = fftphi;
end