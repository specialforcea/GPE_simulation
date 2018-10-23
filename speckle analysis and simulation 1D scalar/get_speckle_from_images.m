
speckle = zeros(49,100,Nx);
M = interpolation_matrix(Nx,Ns);
for i=1:49
    filepath = strcat('speckle bench test data/',num2str(i),'.png');
    
    image = imread(filepath);
    image = double(image(:,1:Ns));
    for j=1:100
        sp = image(500+j,:);
        sp = sp./integr(sp,Ns,deltas);
        inter_sp = (M*sp.').';
        speckle(i,j,:) = inter_sp;
    end
end
    
    