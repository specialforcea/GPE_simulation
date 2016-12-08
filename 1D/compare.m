solve_ODEs;
stoptime = t_ode.*Dip_freq;
length_t = size(t_ode,1);

fp = zeros(length_t,2*order +1);
fp_index = zeros(length_t,2*order +1);

phi_evo = phi_1;

fphi = fourier_transform(phi_evo,Nx);
f_square = (fphi.*conj(fphi))./sum(fphi.*conj(fphi));
[fp(1,:),fp_index(1,:)] = find_peak(f_square,order,k_spacing,Nx);

Deltat = 1e-7;

for i = 2:length_t
    
phi_evo = dynamic(phi_evo,stoptime(i)-stoptime(i-1),Deltat,Beta,Nx,E_r,quanta,k_scale);

fphi = fourier_transform(phi_evo,Nx);
f_square = (fphi.*conj(fphi))./sum(fphi.*conj(fphi));

[fp(i,:),fp_index(i,:)] = find_peak(f_square,order,k_spacing,Nx);
end

pop = cv.*conj(cv);




plot(stoptime,pop(:,order+1),stoptime,f_square(:,order+1))

