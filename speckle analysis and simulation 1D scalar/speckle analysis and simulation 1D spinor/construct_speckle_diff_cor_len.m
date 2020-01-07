function   speckle = construct_speckle_diff_cor_len(image,Nx)
image = double(image);

speckle = zeros(4,30,Nx);
im = image(100:129,:);

for i=1:30
    im(i,:) = im(i,:)-(mean(im(i,:))-std(im(i,:)));
end

im1 = [im,im,im,im,im,im,im];
im1 = im1(:,1:Nx);

speckle(1,:,:) = im1;

im2 = [im,im,im,im];
im2 = im2(:,1:(Nx-1)/2);
M = interpolation_matrix(Nx,(Nx-1)/2);
for i=1:30
    speckle(2,i,:) = (M*im2(i,:).').';
end

im3 = [im,im];
im3 = im3(:,1:(Nx-1)/4);
M = interpolation_matrix(Nx,(Nx-1)/4);
for i=1:30
    speckle(3,i,:) = (M*im3(i,:).').';
end

im4 = im;
im4 = im4(:,1:(Nx-1)/8);
M = interpolation_matrix(Nx,(Nx-1)/8);
for i=1:30
    speckle(4,i,:) = (M*im4(i,:).').';
end