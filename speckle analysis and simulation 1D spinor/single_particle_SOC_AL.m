
OmegaR = 0.5;
band = (-5:0.0001:5);

Deltaq = @(x)(4.*x);
E_up = @(x)(Deltaq(x)/2 + 1/2*sqrt(Deltaq(x).^2+OmegaR^2) + (x-1).^2 );
E_down = @(x)(Deltaq(x)/2 - 1/2*sqrt(Deltaq(x).^2+OmegaR^2) + (x-1).^2 );
dE_up = @(x)(2*x + 8*x./sqrt(Deltaq(x).^2+OmegaR^2));
dE_down = @(x)(2*x - 8*x./sqrt(Deltaq(x).^2+OmegaR^2));

sz = size(band,2);
Ed = E_down(band);
mid = (sz-1)/2;
[~,ind] = min(Ed(mid:end));
min_q = band(mid + ind);
ks = (min_q:0.01:4.2);


a_0 = OmegaR/2;
b_0 = @(x)( Deltaq(x)./2 + 1/2*sqrt(Deltaq(x).^2+OmegaR^2));
c_0 = -a_0;
d_0 = @(x)(-Deltaq(x)./2 + 1/2*sqrt(Deltaq(x).^2+OmegaR^2));

a = @(x)(a_0./sqrt(a_0^2 + b_0(x).^2));
b = @(x)(b_0(x)./sqrt(a_0^2 + b_0(x).^2));
c = @(x)(c_0./sqrt(c_0^2 + d_0(x).^2));
d = @(x)(d_0(x)./sqrt(c_0^2 + d_0(x).^2));

rho = @(x)(2700*( acos(min(1,x./6))- x./6.*sqrt(1-(min(1,x./6)).^2))*2/pi);

scatter_rate_to_up_for = zeros(1,size(ks,2));
scatter_rate_to_up_bac = zeros(1,size(ks,2));
scatter_rate_to_down_for = zeros(1,size(ks,2));
scatter_rate_to_down_bac = zeros(1,size(ks,2));

% % for i=1:size(ks,2)
% %     id_down = find(abs(E_down(band)-E_down(ks(i)))<0.005 & abs(band-ks(i))>0.005 );
% %     
% %     for j=1:size(id_down,2)
% %         scatter_rate_to_down(i) = scatter_rate_to_down(i) + (a(band(id_down(j)))*a(ks(i))+ b(band(id_down(j)))*b(ks(i)))^2*rho(abs(band(id_down(j))-ks(i)));
% %     end
% %     
% %     id_up = find(abs(E_up(band)-E_down(ks(i)))<0.005);
% %     
% %     for j=1:size(id_up,2)
% %         scatter_rate_to_up(i) = scatter_rate_to_up(i) + (c(band(id_up(j)))*a(ks(i))+ d(band(id_up(j)))*b(ks(i)))^2*rho(abs(band(id_up(j))-ks(i)));
% %     end
% %     
% % end

for i=1:size(ks,2)
    
    
    for j=1:size(band,2)
        eng_dif = E_down(band(j))-E_down(ks(i));
        if eng_dif ==0
            fac = 1;
        else
            fac = sin(eng_dif*50)^2/eng_dif^2;
        end
        
        if dE_down(band(j))>=0
            scatter_rate_to_down_for(i) = scatter_rate_to_down_for(i) + fac*(a(band(j))*a(ks(i))+ b(band(j))*b(ks(i)))^2*rho(abs(band(j)-ks(i)));
        else
            scatter_rate_to_down_bac(i) = scatter_rate_to_down_bac(i) + fac*(a(band(j))*a(ks(i))+ b(band(j))*b(ks(i)))^2*rho(abs(band(j)-ks(i)));
            
        end
    end
    
   
    
    for j=1:size(band,2)
        eng_dif = E_up(band(j))-E_down(ks(i));
        if eng_dif ==0
            fac = 1;
        else
            fac = sin(eng_dif*50)^2/eng_dif^2;
        end
        
        if dE_up(band(j))>=0
            scatter_rate_to_up_for(i) = scatter_rate_to_up_for(i) + fac*(c(band(j))*a(ks(i))+ d(band(j))*b(ks(i)))^2*rho(abs(band(j)-ks(i)));
        else
            scatter_rate_to_up_bac(i) = scatter_rate_to_up_bac(i) + fac*(c(band(j))*a(ks(i))+ d(band(j))*b(ks(i)))^2*rho(abs(band(j)-ks(i)));
            
        end
        
    end
    
end

scatter_rate_for = scatter_rate_to_down_for + scatter_rate_to_up_for;
scatter_rate_bac = scatter_rate_to_down_bac + scatter_rate_to_up_bac;
%m = max(scatter_rate);
scatter_rate_for = scatter_rate_for./m1;
scatter_rate_bac = scatter_rate_bac./m2;
% figure(2)


subplot(234)
%plot(band,E_up(band),band,E_down(band))
scatter(band,E_up(band)-min(E_down(band)),[],a(band).^2)

hold on
scatter(band,E_down(band)-min(E_down(band)),[],c(band).^2)
ylim([-0.5,10])
text(0,5,strcat('Omega=',num2str(OmegaR),'E_R'))
ylim([-0.5,10])
ylabel('E/E_R')
xlabel('q/k_R')

subplot(235)

%scatter(scatter_rate,E_down(ks))
hold on
plot(scatter_rate_bac,E_down(ks)-min(E_down(band)))
ylim([-0.5,10])
xlabel('normalized scatter rate')

subplot(236)

%scatter(scatter_rate,E_down(ks))
hold on
plot(scatter_rate_for,E_down(ks)-min(E_down(band)))
ylim([-0.5,10])
xlabel('normalized scatter rate')



%single spin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = (0.0:0.01:3.2);
band = (-5:0.0001:5);
scatter_rate_single_for = zeros(1,size(ks,2));
scatter_rate_single_bac = zeros(1,size(ks,2));
E0 = @(x)(x.^2);
for i=1:size(ks,2)
    
    
    for j=1:size(band,2)
        eng_dif = E0(band(j))-E0(ks(i));
        if eng_dif ==0
            fac = 1;
        else
            fac = sin(eng_dif*50)^2/eng_dif^2;
        end
        if j<size(band,2)/2
            
            scatter_rate_single_for(i) = scatter_rate_single_for(i) + fac*rho(abs(band(j)-ks(i)));
        else
            scatter_rate_single_bac(i) = scatter_rate_single_bac(i) + fac*rho(abs(band(j)-ks(i)));
        end
    end
    
end

m1 = max(scatter_rate_single_for);
scatter_rate_single_for = scatter_rate_single_for/m1;
m2 = max(scatter_rate_single_bac);
scatter_rate_single_bac = scatter_rate_single_bac/m2;


subplot(231)

scatter(band,E0(band),10,'o')
ylim([-0.5,10])
text(0,5,strcat('Omega=0'))
ylabel('E/E_R')
xlabel('q/k_R')
subplot(232)

%scatter(scatter_rate_single,E0(ks))
hold on
plot(scatter_rate_single_for,E0(ks))
xlabel('normalized scatter rate')
ylim([-0.5,10])

subplot(233)

%scatter(scatter_rate_single,E0(ks))
hold on
plot(scatter_rate_single_bac,E0(ks))
xlabel('normalized scatter rate')
ylim([-0.5,10])