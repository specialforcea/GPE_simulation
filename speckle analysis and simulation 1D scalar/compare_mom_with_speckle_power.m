
fspeckle = zeros(5,Nx);
for i=8:12
    speckle = load(strcat('speckle bench test data/12',num2str(i),'.mat'));
    speckle = speckle.inten/1e5;
    speckle = speckle - repmat(mean(speckle,2),1,Nx);
    
    fsp = fourier_transform(speckle',Nx,deltax)';
    fsp = abs(fsp);
    fsp = mean(fsp);
    fspeckle(i-7,:) = fsp;
end
    
    
    