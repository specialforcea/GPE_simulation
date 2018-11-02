li=[20,25,30];
fspeckle = zeros(3,Ns);
for i=1:3
    speckle = imread(strcat('speckle bench test data/12/',num2str(li(i)),'mm.png'));
    speckle = double(speckle);
    speckle = speckle - repmat(mean(speckle,2),1,Ns);
    
    fsp = fourier_transform(speckle',Ns,deltas)';
    sq_fsp = sq(fsp);
    
    power_sprectra = abs(inverse_ft(sq_fsp',Ns,deltasf)');
    
    pow = mean(power_sprectra);
    
    fspeckle(i,:) = pow;
end
    
    
    