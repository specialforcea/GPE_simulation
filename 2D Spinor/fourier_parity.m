function [paritx,parity] = fourier_parity(Nx,Ny)

lin_x = linspace(0,Nx-1,Nx);
unitx = ones(Ny,1);
lin_y = linspace(0,Ny-1,Ny);
unity = ones(1,Nx);


% lin_x = linspace(0,Nx-2,Nx-1);
% unitx = ones(Ny-1,1);
% lin_y = linspace(0,Ny-2,Ny-1);
% unity = ones(1,Nx-1);

px = (-1).^lin_x;
px2 = unitx*px;
paritx = px2;
paritx(:,:,2) = px2;
paritx(:,:,3) = px2;


py = (-1).^lin_y;
py2 = py'*unity;
parity = py2;
parity(:,:,2) = py2;
parity(:,:,3) = py2;

