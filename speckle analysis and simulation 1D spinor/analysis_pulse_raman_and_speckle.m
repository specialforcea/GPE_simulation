file_path = '/Users/yueyuchen/Documents/Academy/Research/GPE/yuchen_GPE/GPE_simulation/speckle analysis and simulation 1D spinor/simulation_results/12082020pulse_raman_and_speckle';


load(strcat(file_path,'/6/final_phi.mat'));


fphi1 = zeros(2,Nx);


i = 10;

fp = reshape(final_phi(i,:,:,:),[20,2,Nx]);
for j=1:20
fp1 = reshape(fp(j,:,:),[2,Nx]);
fp1 = fourier_transform(fp1,Nx,deltax);
fp1 = sq(fp1);
fphi1 = fphi1 + fp1;
end

fphi1 = fphi1./20;
plot(f/k_R*2*pi, fphi1(1,:))
hold on

xlabel('k/k_R')
ylabel('momentum distribution')
% legend('\Omega=1','\Omega=3','\Omega=6')