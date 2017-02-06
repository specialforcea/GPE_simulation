function dispersion = fourier_dispersion(Nx,Ny,Lx,Ly)
%phi_up = phi_up/norm1d(phi_up, Nx, DeltaX);
% if(mod(Nx,2) ~= 0)
%     lin_n =  (0:1:Nx-1) - (Nx-1)/2;
% else
    %N = Nx + 1;
    lin_nx = (0:1:Nx-1) - (Nx-1)/2;
    lin_nx = lin_nx(1:Nx-1);
%end
% if(mod(Nx,2) ~= 0)
%     lin_n =  (0:1:Nx-1) - (Nx-1)/2;
% else
%     N = Nx + 1;
    lin_ny = (0:1:Ny-1) - (Ny-1)/2;
    lin_ny = lin_ny(1:Ny-1);
% end
dispx = -(2*pi/Lx.*lin_nx).^2./2;
unity = ones(Ny-1,1);
disp2 = unity*dispx;
dispersionx = disp2;
dispersionx(:,:,2) = disp2;
dispersionx(:,:,3) = disp2;

dispy = -(2*pi/Ly.*lin_ny).^2./2;
unity = ones(1,Nx-1);
disp2 = dispy'*unity;
dispersiony = disp2;
dispersiony(:,:,2) = disp2;
dispersiony(:,:,3) = disp2;

dispersion = dispersionx + dispersiony;
