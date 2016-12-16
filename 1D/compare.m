solve_ODEs;
stoptime = t_ode.*Dip_freq;
length_t = size(t_ode,1);

fp = zeros(length_t,2*order +1);
fp_index = zeros(length_t,2*order +1);

phi_evo = phi_0;

fphi = fourier_transform(phi_evo,Nx,deltax);

[fp(1,:),fp_index(1,:)] = find_peak(sq(fphi),order,k_spacing,Nx);



Deltat = 1e-7;   

for i = 2:length_t
 
phi_evo = dynamic(phi_evo,stoptime(i)-stoptime(i-1),Deltat,Beta,Nx,V,k_scale,f,deltax,deltaf,L);

fphi = fourier_transform(phi_evo,Nx,deltax);


[fp(i,:),fp_index(i,:)] = find_peak(sq(fphi),order,k_spacing,Nx);
end
s = sum(fp,2);
div = zeros(size(fp));
for i = 1:2*order+1
    div(:,i) = s;
end
fp = fp./div;
pop = cv.*conj(cv);




%plot(stoptime,pop(:,order+1),stoptime,f_square(:,order+1))

