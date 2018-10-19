Deltat = 1e-9;
Stop_crit = 8.47e-10;



potential = @(x,y)(1/2*x.^2 + 1/2*gamma^2*y.^2);


init_spin = [1 0 0];
phi_0 = Thomas_fermi2D(c0,X,Y,TF_radius,potential(X,Y),gamma);
phi = phi_0;

fftNx = Nx -1;
fftNy = Ny -1;
fftphi = phi(1:fftNx,1:fftNy);
fftX = X(1:fftNx,1:fftNy);
fftY = Y(1:fftNx,1:fftNy);
fftL = Lx - deltax;
fftfx = fx(1:fftNx);

difference = 1;
evo = 200;
n = 0;
while (difference)
    
    phi_up = imtime_evolve(fftphi, potential, Deltat, V0,k_scale,fftX,fftY, fftNx,fftNy, deltax,deltay,deltafx,deltafy,c0,paritx(1:fftNx,1:fftNy,:),parity(1:fftNx,1:fftNy,:),dispersion);
    if (max(abs(phi_up(:)-fftphi(:))) < Stop_crit)
        difference = 0;
    end
    
    if (mod(n,evo) == 0)
        max(abs(phi_up(:)-fftphi(:)))
%         
    end
    fftphi = phi_up;
    n = n + 1;
    
end
phi_up = zeros(Nx,Ny);
phi_up(1:fftNx,1:fftNy) = fftphi;
phi_up(Nx,:) = phi_up(1,:);
phi_up(:,Ny) = phi_up(:,1);
% r1 = randn(Nx,Ny) + 1i*randn(Nx,Ny);
% r1 = r1/norm2d(r1,Nx,Ny,deltax,deltay);
% r2 = r1 - integr2d(conj(phi_up).*r1,Nx,Ny,deltax,deltay).*phi_up;
% phi_up = phi_up + 0.3*r2;



phi_0(:,:,2) = init_spin(2).*phi_up;
phi_0(:,:,3) = init_spin(3).*phi_up;
phi_0(:,:,1) = init_spin(1).*phi_up;

phi_0 = phi_0./norm2d(phi_0,Nx,Ny,deltax,deltay);




