filepath = strcat('speckle bench test data/',num2str(1),'.png');
image = imread(filepath);
image = double(image);
speckle = construct_speckle_diff_cor_len(image,Nx);

phi_all = zeros(4,6,Nx);
intensity = [0.25,0.33,0.5,1,2,3];
for i=1:4
    for j=1:6
        phi_tem = zeros(1,Nx);
        sp = speckle(i,:,:).*intensity(j);
        sp = reshape(sp,[30,Nx]);
        for h=1:30
            phi = dynamic(phi_0,12.0,1e-3,c0,c2,Nx,sp(h,:),k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);
            phi_tem = phi_tem + sq(phi(1,:));
        end
        phi_all(i,j,:) = phi_tem./30;
    end
end

