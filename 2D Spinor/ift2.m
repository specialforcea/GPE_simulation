function y = ift2(phi,paritx,parity,deltafx,deltafy,Nx,Ny)



y = zeros(Ny,Nx,3);

for i = 1:Ny
    y(i,:,1) = inverse_ft(phi(i,:,1),Nx,deltafx);
    y(i,:,2) = inverse_ft(phi(i,:,2),Nx,deltafx);
    y(i,:,3) = inverse_ft(phi(i,:,3),Nx,deltafx);
end

for i = 1:Nx
    y(:,i,1) = inverse_ft(y(:,i,1)',Ny,deltafy)';
    y(:,i,2) = inverse_ft(y(:,i,2)',Ny,deltafy)';
    y(:,i,3) = inverse_ft(y(:,i,3)',Ny,deltafy)';
end


end