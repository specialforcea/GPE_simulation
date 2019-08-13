%hbar^2k_R^2/2/m~2pi*3.7kHz, scale it to 3.7
%mean(v(x))=2pi*1kHz, scale it to 1
i=8;
filepath = strcat('speckle bench test data/numerical_speckle/13/inten_', num2str(i),'.mat');
speckle = load(filepath);
speckle = speckle.inten;

speckle = speckle/mean(mean(speckle));%average intensity about 1 in simulation units

speckle = speckle - mean(mean(speckle));

speckle_k = mean(fourier_transform(speckle,Nx,deltax));


rev = fliplr(speckle_k)';
%plot(f./k_R*2*pi,sq(speckle_k))

kine = diag((f./k_R*2*pi).^2*3.7);

V = zeros(Nx,Nx);

half = (Nx-1)/2;

V(:,half) = rev;


for i=1:half
    V(1:Nx-i,half-i) = rev(1+i:end);
    V(1+i:end,half+i) = rev(1:Nx-i);
end





