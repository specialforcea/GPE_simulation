function phi_up = time_evolve(phi, potential, speckle,Deltat, X, Nx, deltax,deltaf,L,c0,c2,Omega,k_R,detuning)

density = sq(phi);
tot_density = sum(density,1);

p = potential(X) + speckle;
pot = [p;p;p];

nonlin_pot = zeros(3,Nx);
nonlin_pot(1,:) = c0.*tot_density + c2.*(density(1,:) + density(2,:) - density(3,:));
nonlin_pot(2,:) = c0.*tot_density + c2.*(density(1,:) + density(3,:)) + 0*detuning;
%nonlin_pot(3,:) = c0.*tot_density + c2.*(density(3,:) + density(2,:) - density(1,:)) - 500*detuning;



diag_pot = -Deltat.*(pot + nonlin_pot);

exp_pot = exp(diag_pot);

phi_up = exp_pot.*phi;

if(mod(Nx,2) ~= 0)
    lin_n =  (0:1:Nx-1) - (Nx-1)/2;
else
    N = Nx + 1;
    lin_n = (0:1:N-1) - (N-1)/2;
    lin_n = lin_n(1:Nx);
end
dispersion = -(2*pi/L*lin_n).^2 *Deltat/2;
evol = exp(dispersion);

fourier_phi = zeros(3,Nx);
fourier_phi(1,:) = fourier_transform(phi_up(1,:),Nx,deltax);
fourier_phi(2,:) = fourier_transform(phi_up(2,:),Nx,deltax);
%fourier_phi(3,:) = fourier_transform(phi_up(3,:),Nx,deltax);

fourier_phi_evo = [evol.*fourier_phi(1,:);evol.*fourier_phi(2,:);evol.*fourier_phi(3,:)];

phi_up(1,:) = inverse_ft(fourier_phi_evo(1,:),Nx,deltaf);
phi_up(2,:) = inverse_ft(fourier_phi_evo(2,:),Nx,deltaf);
%phi_up(3,:) = inverse_ft(fourier_phi_evo(3,:),Nx,deltaf);

phi_up = phi_up/norm1d(phi_up, Nx, deltax);

CA = -Omega*exp(-1i*2*k_R*X) + c2*phi_up(2,:).*0.*conj(phi_up(3,:));%coupling between 1 and 0
CC = -Omega*exp(-1i*2*k_R*X) + c2*phi_up(1,:).*conj(phi_up(2,:));%coupling between -1 and 0

% CA = -Omega + c2*phi_up(2,:).*0.*conj(phi_up(3,:));%coupling between 1 and 0
% CC = -Omega + c2*phi_up(1,:).*conj(phi_up(2,:));%coupling between -1 and 0


MatV = CC./((sqrt(1+sq(CC)./sq(CA))).*conj(CA));
MatV(2,:) = zeros(1,Nx);
MatV(3,:) = 1./sqrt(1 + sq(CC)./sq(CA));
MatV(1,:,2) = CA./(sqrt(2*(sq(CA) + sq(CC))./sq(CC)).*conj(CC));
MatV(2,:,2) = -sqrt(sq(CC))./(sqrt(2).*conj(CC));
MatV(3,:,2) = sqrt(sq(CC))./sqrt(2.*(sq(CA)+sq(CC)));
MatV(1,:,3) = CA./(sqrt(2*(sq(CA) + sq(CC))./sq(CC)).*conj(CC));
MatV(2,:,3) = sqrt(sq(CC))./(sqrt(2).*conj(CC));
MatV(3,:,3) = sqrt(sq(CC))./sqrt(2.*(sq(CA)+sq(CC)));

MatV = permute(MatV,[1,3,2]);
MatB = permute(MatV,[2,1,3]);
MatB = conj(MatB);

eigens = zeros(3,Nx);
eigens(1,:) = zeros(1,Nx);
eigens(2,:) = -sqrt(sq(CA) + sq(CC));
eigens(3,:) = sqrt(sq(CA) + sq(CC));

% for i = 1:Nx
%     phi_up(:,i) = MatV(:,:,i)*(exp(-eigens(:,i).*Deltat).*(MatB(:,:,i)*phi_up(:,i)));
% end

temp = zeros(3,Nx);

for i=1:3
	temp(i,:) = sum(reshape(MatB(i,:,:),[3,Nx]).*phi_up).*exp(-eigens(i,:).*Deltat);
end

for i=1:2
	phi_up(i,:) = sum(reshape(MatV(i,:,:),[3,Nx]).*temp);
end




phi_up = phi_up/norm1d(phi_up, Nx, deltax);

end



