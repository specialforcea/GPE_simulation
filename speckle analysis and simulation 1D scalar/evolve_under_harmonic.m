mk = 1;
phi_mk = phi_0;
phi_mk(1,:) = phi_mk(1,:).*exp(1i*k_R*mk*X);


phi = phi_mk;

i=8;
filepath = strcat('speckle bench test data/numerical_speckle/13/inten_', num2str(i),'.mat');
speckle = load(filepath);
speckle = speckle.inten;

speckle = speckle/1e6;%average intensity about 1 in simulation units

speckle = speckle*5;%make it 68% of chemical potential.


for i=1:700
    phi_1 = dynamic(phi,0.01,1e-4,0,0,Nx,0,1,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);

    fp = sq(fourier_transform(phi_1(1,:),Nx,deltax));
    
    plot(f./k_R*2*pi,fp);
    title(num2str(i*0.01/Dip_freq));
    drawnow;

    phi = phi_1;

end


