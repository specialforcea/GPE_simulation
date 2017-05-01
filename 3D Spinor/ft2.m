function y = ft2(phi,paritx,parity,deltax,deltay)

Nx = size(phi,2);
Ny = size(phi,1);

y = zeros(Ny,Nx,3);
for i = 1:Ny
    y(i,:,1) = fourier_transform(phi(i,:,1),Nx,deltax);
    y(i,:,2) = fourier_transform(phi(i,:,2),Nx,deltax);
    y(i,:,3) = fourier_transform(phi(i,:,3),Nx,deltax);
end

for i = 1:Nx
    y(:,i,1) = fourier_transform(y(:,i,1)',Ny,deltay)';
    y(:,i,2) = fourier_transform(y(:,i,2)',Ny,deltay)';
    y(:,i,3) = fourier_transform(y(:,i,3)',Ny,deltay)';
end

end