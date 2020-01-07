
phi = phi_0;
popu = zeros(3,100);

for i=1:100
    phi_1 = dynamic(phi,1e-3,1e-5,c0,c2,Nx,0,0,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax,k_R,detuning);
    popu(1,i) = integr(sq(phi_1(1,:)),Nx,deltax);
    popu(2,i) = integr(sq(phi_1(2,:)),Nx,deltax);
    popu(3,i) = integr(sq(phi_1(3,:)),Nx,deltax);
    
    phi = phi_1;
%     plot(X,sq(phi_1(1,:)))
%     drawnow
end