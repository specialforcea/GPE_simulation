function f = ODEs_define(t,c)
order = 1;%to nth order
hbar = 1.05*10^-34;
Rb_Mass = 1.42*10^-25;% mass of Rb87 in Kg
raman_wavelength = 790*10^-9; % raman beam wave length
k = 2*pi/raman_wavelength; %raman beam wave vector
E_r = (hbar)^2*(2*k)^2/(2*Rb_Mass);%recoil energy
lattice_pot = 1*E_r;
A = -1i*E_r/hbar;
B = -1i*lattice_pot/(4*hbar);
dcdt(1) = (A*order^2+2*B)*c(1)+B*c(2);
dcdt(2*order+1) = B*c(2*order)+(A*order^2+2*B)*c(2*order+1);

for i = 2:2*order

dcdt(i) = B*c(i-1)+(A*(order-i+1)^2+2*B)*c(i)+B*c(i+1);

end
f = zeros(2*order+1,1);
for i = 1:2*order+1
    f(i) = dcdt(i);

end