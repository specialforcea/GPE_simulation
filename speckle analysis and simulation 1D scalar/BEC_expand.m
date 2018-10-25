mom_evo = zeros(5,300);

phi = phi_0;

for i=1:300
    phi_1 = dynamic(phi,1e-2,1e-5,c0,c2,Nx,0,0,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);

    fp = sq(fourier_transform(phi_1(1,:),Nx,deltax));
    


    mom_evo(1,i) = sqrt(integr(f.^2.*fp,Nx,deltaf));

            plot(f(1,3500:4500)./k_spacing,fp(1,3500:4500))
            drawnow

    phi = phi_1;

end