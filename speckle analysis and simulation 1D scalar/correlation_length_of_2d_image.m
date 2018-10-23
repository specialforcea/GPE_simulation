function [xcorr,ycorr] = correlation_length_of_2d_image(image)
[Ny,Nx] = size(image);
xcorr = zeros(1,Nx);
ycorr = zeros(1,Ny);

avg = mean(mean(image));
image = image - avg;

Var = sum(sum(image.*image))/Nx/Ny;
xcorr(1) = 1;
ycorr(1) = 1;
for i=2:Nx
    tem_corr = 0;
    j = 1;
    while j+i-1<=Nx
        tem_corr = tem_corr + image(:,j).*image(:,j+i-1);
        j = j + 1;
    end
    xcorr(i) = sum(tem_corr)/((Nx-i+1)*Ny-1)/Var;
end

for i=2:Ny
    tem_corr = 0;
    j = 1;
    while j+i-1<=Ny
        tem_corr = tem_corr + image(j,:).*image(j+i-1,:);
        j = j + 1;
    end
    ycorr(i) = sum(tem_corr)/((Ny-i+1)*Nx-1)/Var;
end
end