function phi_2 = dynamic(phi,t_evo,Deltat,Beta,Nx,E_r,quanta,k_scale,f)

 
Stop_time = t_evo;


TF_radius = (3*Beta/2)^(1/3);
xmin = -TF_radius;
xmax = TF_radius;



lattice_pot = 4*E_r;
V = lattice_pot/quanta;
DeltaX = (xmax-xmin)/(Nx-1);


potential = @(x)(1/2*x.^2 + V.*sin(k_scale.*x).^2 );
%order = 2;

X = linspace(xmin,xmax,Nx);
%N_tf = int16(TF_radius/DeltaX);
%cutX = X < TF_radius;

t = 0;
fftNx = Nx -1;
fftphi = phi(1:fftNx);
fftX = X(1:fftNx);
Deltat = 1i*Deltat;
evo = 500;
n = 0;
while (t < Stop_time)
    fftphi = time_evolve(fftphi, potential,Deltat,fftX,Beta,fftNx,DeltaX);
    t = t - 1i*Deltat;
    n = n + 1;
    
    if (mod(n,evo)==0)
        fphi = fourier_transform(fftphi,fftNx);
        b = fphi.*conj(fphi);
        plot(f(1:fftNx),b);
        drawnow;
        
    end
    
end
phi_2 = zeros(1,Nx);
phi_2(1:fftNx) = fftphi;
phi_2(Nx) = phi_2(1);

end