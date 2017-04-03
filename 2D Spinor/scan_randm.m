


randm = [0.1 0.3 0.9 1];
ratio = 1;
%phibig = zeros(Nx,Ny,12);
for j = 1:4
phi_up = zeros(Nx,Ny);
phi_up(1:fftNx,1:fftNy) = fftphi;
phi_up(Nx,:) = phi_up(1,:);
phi_up(:,Ny) = phi_up(:,1);
r1 = randn(Nx,Ny) + 1i*randn(Nx,Ny);
r1 = r1/norm2d(r1,Nx,Ny,deltax,deltay);
r2 = r1 - integr2d(conj(phi_up).*r1,Nx,Ny,deltax,deltay).*phi_up;
phi_up = phi_up + randm(j)*r2;



phi_0(:,:,2) = init_spin(2).*phi_up;
phi_0(:,:,3) = init_spin(3).*phi_up;
phi_0(:,:,1) = init_spin(1).*phi_up;

phi_0 = phi_0./norm2d(phi_0,Nx,Ny,deltax,deltay);

pop_evolve;
%phibig(:,:,j) = phi_evo(:,:,1);
% ff = fourier_transform2(phi_evo,paritx,parity,deltax,deltay);
% ratio = integr2d(sq(ff(128-4:128+4,128-4:128+4,1)),9,9,deltafx,deltafy);
% avphi = mean(phi_evo((Nx-1)/2-10:(Nx-1)/2+10,:,1));
% avphi0 = mean(phi_0((Nx-1)/2-10:(Nx-1)/2+10,:,1));
% figure;
% plot(X(1,:),sq(avphi),X(1,:),sq(avphi0));
% title(randm);
% savefig(strcat(num2str(randm),'.fig'));
% if ratio >  0.9
%     randm = randm + 0.001;
% end
end



