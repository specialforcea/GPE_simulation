figure(5)
for O=1:4
T = 1;

phi = reshape(final_phi(O,T,:,:,:),[20,2,Nx]);
fphi = zeros(2,Nx);
for i=1:20
    ffphi = fourier_transform(reshape(phi(i,:,:),[2,Nx]),Nx,deltax);
    fphi = fphi + sq(ffphi);
end
fphi = fphi./20;

plot(f./k_R*2*pi,fphi(1,:))
ylim([0,0.001])
hold on
end



% T = 1;
% 
% phi = reshape(final_phi(T,:,:,:),[20,2,Nx]);
% fphi = zeros(2,Nx);
% for i=1:20
%     ffphi = fourier_transform(reshape(phi(i,:,:),[2,Nx]),Nx,deltax);
%     fphi = fphi + sq(ffphi);
% end
% fphi = fphi./20;
% 
% plot(f./k_R*2*pi,fphi(1,:))
% ylim([0,0.001])

