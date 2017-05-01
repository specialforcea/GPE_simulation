function dispersion = fourier_dispersion(Nx,Ny,Nz,Lx,Ly,Lz)
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


lin_nz = (0:1:Nz-1) - (Nz-1)/2;
lin_nz = lin_nz(1:Nz-1);

dispx = -(2*pi/Lx.*lin_nx).^2./2;
unity = ones(Ny-1,1);
disp1 = unity*dispx;
% dispersionx = disp2;
% dispersionx(:,:,2) = disp2;
% dispersionx(:,:,3) = disp2;

dispy = -(2*pi/Ly.*lin_ny).^2./2;
unity = ones(1,Nx-1);
disp2 = dispy'*unity;
% dispersiony = disp2;
% dispersiony(:,:,2) = disp2;
% dispersiony(:,:,3) = disp2;
dispz = -(2*pi/Lz.*lin_ny).^2./2;

disp3 = disp1 + disp2;
disp = zeros(Ny-1,Nx-1,Nz-1);
for i = 1:Nz-1
	disp(:,:,i) = disp3 + dispz(i);
end
dispersion = disp;
dispersion(:,:,:,2) = disp;
dispersion(:,:,:,3) = disp;
