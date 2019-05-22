phi_ei = zeros(3,Nx);
phi_ei(1,:) = phi_0(1,:).*cost(1);
phi_ei(2,:)=  phi_0(1,:).*sint(1).*exp(1i.*2.*k_R.*X);
phi = phi_ei.*[exp(1i.*kick.*k_R.*X);exp(1i.*kick.*k_R.*X);exp(1i.*kick.*k_R.*X)];
spin1 = zeros(1,320);
draw = 0;
    
for i=1:320
    phi_dress = dynamic(phi,0.005,1e-5,c0,c2,Nx,0,0,0,k_scale,f,deltax,deltaf,L,OmegaR*detuning,xmin,xmax,k_R,0);
    phi = phi_dress;
    spin1(i) = integr(sq(phi(1,:)),Nx,deltax);

    if draw==1 && mod(i,5)==0
        subplot(121)
        plot(1:1:320,spin1)
        subplot(122)
        plot(X,sq(phi))
        drawnow;
    end
    
end