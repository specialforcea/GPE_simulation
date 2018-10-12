init_phi = sqrt(gaussian(X,0,0.04));%500nK
init_phi = [init_phi;0*init_phi;0*init_phi];
phi_2 = dynamic(init_phi,13.0,1e-5,c0,c2,Nx,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);
