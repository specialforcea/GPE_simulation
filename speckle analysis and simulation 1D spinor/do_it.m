ks = (1.2:0.1:3.2);
OmegaR = 0.5;
band = (-4:0.001:4);
y1 = (band+1).^2;
y2 = (band-1).^2;
dt = abs(y2-y1);
up = min(y1,y2) + 1/2*(dt+sqrt(dt.^2+OmegaR^2));
down = min(y1,y2) + 1/2*(dt-sqrt(dt.^2+OmegaR^2));

k_minima = sqrt((4*ks).^2./(OmegaR^2 + (4*ks).^2));
dtn = 4.*(ks - k_minima);

savepath = 'simulation_results/05012019SOC_dressed_state_8TF13/phi_1_50ms_';

for j=1:21
save(strcat(savepath,num2str(ks(j)),'k_R_',num2str(OmegaR),'E.mat'),'phi_1');

kick = ks(j) - k_minima(j);
phi = phi_1.*[exp(1i.*kick.*k_R.*X);exp(1i.*kick.*k_R.*X);exp(1i.*kick.*k_R.*X)];
end