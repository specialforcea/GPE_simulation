h = 6.63e-34;
x = (-250:0.01:250)*1e-6;
P1 = 1;
P2 = (0.01:0.01:0.8);
w1 = 230e-6;
w2 = 200e-6;
u = 50e-6;
I1 = 2*P1/pi/w1^2.*exp(-2*x.^2./w1^2);
I2 = 2/w2^2.*exp(-2*(x-u).^2./w2^2);
V = -h*0.003*I1;
E_L = h*300;%Hz,at raman angel 33.
minima = zeros(1,size(P2,2));
n = (0.1:0.1:2.5);
xn = sqrt(1/2*(-log(1-n.^2./120))).*w1;
figure(1)
plot(n,xn)
for i=1:size(P2,2)
I = -I1-P2(i).*I2;
[a,ind] = min(I);
minima(i) = x(ind);
end

figure(2)
plot(P2,minima)


    