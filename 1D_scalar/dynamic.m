function phi_2= dynamic(phi,t_evo,Deltat,c0,Nx,trap,V,k_scale,f,deltax,deltaf,L,xmin,xmax)

 
Stop_time = t_evo;










potential = @(x)(trap/2*x.^2 + V.*sin(k_scale.*x).^2 );


X = linspace(xmin,xmax,Nx);
t = 0;
fftNx = Nx -1;
fftphi = phi(1:fftNx);
fftX = X(1:fftNx);
fftL = L - deltax;
fftf = f(1:fftNx);


Deltat = 1i*Deltat;

%show some results every evo of evolutions
evo = 5000;
n = 0;
draw = 1;

m=1;
while (t < Stop_time)
                          
    fftphi = time_evolve(fftphi, potential,Deltat,fftX,fftNx,deltax,deltaf,fftL,c0);
    t = t - 1i*Deltat;
    n = n + 1;
   
    
    if (mod(n,evo)==0 && draw == 1)       
        m = m + 1;
        plot(fftX,sq(fftphi))
        
        drawnow;
        
    end
    
end
phi_2 = zeros(1,Nx);
phi_2(1:fftNx) = fftphi;
phi_2(Nx) = phi_2(1);


end