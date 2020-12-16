file_path = '/Users/yueyuchen/Documents/Academy/Research/GPE/yuchen_GPE/GPE_simulation/speckle analysis and simulation 1D spinor/simulation_results/09172020SOC_speckle_pulse_200ms_cluster';


load(strcat(file_path,'/final_phi.mat'));


fphi1 = zeros(4,2,Nx);


i = 10;

for o=1:4
fp = reshape(final_phi(o,i,:,:,:),[20,2,Nx]);

for j=1:20
fp1 = reshape(fp(j,:,:),[2,Nx]);
fp1 = fourier_transform(fp1,Nx,deltax);
fp1 = reshape(sq(fp1),[1,2,Nx]);
fphi1(o,:,:) = fphi1(o,:,:) + fp1;
end
end
fphi1 = fphi1./20;

for o=1:4
plot(f/k_R*2*pi, reshape(fphi1(o,2,:),[1,Nx]))
hold on
end

xlabel('k/k_R')
ylabel('momentum distribution')
% legend('\Omega=1','\Omega=3','\Omega=6')
% fphi(1:4,:,:) = fphi(1:4,:,:) + fphi1;