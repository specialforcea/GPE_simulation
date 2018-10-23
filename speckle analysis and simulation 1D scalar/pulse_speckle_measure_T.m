speckle = load('speckle bench test data\T_simu\inten_12.mat');
speckle = speckle.inten;
speckle = speckle./1e5;%average intensity about 5 in simulation units
speckle = speckle(31:60,:);
v_width = zeros(6,30);
chem1 = zeros(6,30);
T1 = zeros(6,30);
speckle_pulse_time = 0.88;%2ms
speckle_factor = [1, 2, 5, 10, 16.67, 26.67];
for i=1:6
    speckle1 = speckle.*speckle_factor(i);
    for j=1:30
    
    speckle_row = speckle1(j,:);
    trap = 1;%dipole trap on
    phi_1 = dynamic(phi_0,speckle_pulse_time,1e-4,c0,c2,Nx,speckle_row,trap,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);
    
    hold_time = 2.2;%5ms
    phi_2 = dynamic(phi_1,hold_time,1e-4,c0,c2,Nx,0,trap,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);
    
    
    
%     trap = 0;
%     tof = 11;%25ms
%     phi_3 = dynamic(phi_2,tof,1e-4,c0,c2,Nx,0,trap,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);
    fphi = fourier_transform(phi_2(1,:),Nx,deltax);
    v_width(i,j) = sqrt(integr(f.^2.*sq(fphi),Nx,deltaf));
    chem1(i,j) = chem_pot(phi_2(1,:),X,Nx,c0,k_scale,deltax,deltaf,V,L);
    T1(i,j)  = (v_width(i,j)*a_0*Dip_freq)^2*Rb_Mass/k_b;
    
    end
end