function [eigens,MatB,MatV] = compute_omega_matrix(Omega,k_R,X,Nx)
CA = -Omega/2*exp(-1i*2*k_R*X); %+ c2*phi_up(2,:).*conj(phi_up(3,:));%coupling between 1 and 0
CC = -Omega/2*exp(-1i*2*k_R*X); %+ c2*phi_up(1,:).*conj(phi_up(2,:));%coupling between -1 and 0

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
end