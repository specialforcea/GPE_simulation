[Ny,Nx] = size(image);
[parity_x,parity_y] = fourier_parity(Nx,Ny);
parity_x = parity_x(:,:,1);
parity_y = parity_y(:,:,1);
fimage = fourier_transform2(image,parity_x,parity_y,1,1);
cimage = inverse_ft2(sq(fimage),parity_x,parity_y,1,1,Nx,Ny);
imshow(cimage)