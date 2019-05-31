ks = (1.2:0.01:5.0);
OmegaR = 4.0;
band = (-5:0.0001:5);

Deltaq = @(x)(4.*x);
E_up = @(x)(Deltaq(x)/2 + 1/2*sqrt(Deltaq(x).^2+OmegaR^2) + (x-1).^2 );
E_down = @(x)(Deltaq(x)/2 - 1/2*sqrt(Deltaq(x).^2+OmegaR^2) + (x-1).^2 );



a_0 = OmegaR/2;
b_0 = @(x)( Deltaq(x)./2 + 1/2*sqrt(Deltaq(x).^2+OmegaR^2));
c_0 = -a_0;
d_0 = @(x)(-Deltaq(x)./2 + 1/2*sqrt(Deltaq(x).^2+OmegaR^2));

a = @(x)(a_0./sqrt(a_0^2 + b_0(x).^2));
b = @(x)(b_0(x)./sqrt(a_0^2 + b_0(x).^2));
c = @(x)(c_0./sqrt(c_0^2 + d_0(x).^2));
d = @(x)(d_0(x)./sqrt(c_0^2 + d_0(x).^2));

rho = @(x)(( acos(min(1,x./6))- x./6.*sqrt(1-(min(1,x./6)).^2))*2/pi);

scatter_rate_to_up = zeros(1,size(ks,2));
scatter_rate_to_down = zeros(1,size(ks,2));
for i=1:size(ks,2)
    id_down = find(abs(E_down(band)-E_down(ks(i)))<0.005 & band ~= ks(i));
    
    for j=1:size(id_down,2)
        scatter_rate_to_down(i) = scatter_rate_to_down(i) + (a(band(id_down(j)))*a(ks(i))+ b(band(id_down(j)))*b(ks(i)))^2*rho(abs(band(id_down(j))-ks(i)));
    end
    
    id_up = find(abs(E_up(band)-E_down(ks(i)))<0.005);
    
    for j=1:size(id_up,2)
        scatter_rate_to_up(i) = scatter_rate_to_up(i) + (c(band(id_up(j)))*a(ks(i))+ d(band(id_up(j)))*b(ks(i)))^2*rho(abs(band(id_up(j))-ks(i)));
    end
    
end

scatter_rate = scatter_rate_to_down + scatter_rate_to_up;
scatter_rate = scatter_rate./max(scatter_rate);

subplot(121)
scatter(ks,scatter_rate)
hold on
plot(ks,scatter_rate)
ylabel('normalized scatter rate')
xlabel('q_0/k_R')
subplot(122)
plot(band,E_up(band),band,E_down(band))
text(0,5,strcat('Omega=',num2str(OmegaR),'E_R'))
ylabel('E/E_R')
xlabel('q/k_R')