scale_parameters;
Groundstate;

i=8;
filepath = strcat('speckle bench test data/numerical_speckle/13/inten_', num2str(i),'.mat');
speckle = load(filepath);
speckle = speckle.inten;

speckle = speckle/1e6;%average intensity about 1 in simulation units

speckle = speckle*5;%make it 68% of chemical potential.

save_path = 'simulation_results/09052019pulse_speckle/';

intensity = [1,2,5,10,15,20];

%momentum_dist1 = zeros(6,20,25,Nx);
wid = zeros(6,20,125);

for i=1:size(intensity,2)
    rand_row = randi([1 Nx],1,20);
    for j=1:20

        speckle_row = speckle(rand_row(j),:).*intensity(i);

        phi = phi_0;

        for k=1:125
            phi_1 = dynamic(phi,0.0005,1e-5,c0,c2,Nx,speckle_row,1,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);
            
            
                
            fp = sq(fourier_transform(phi_1(1,:),Nx,deltax));
            %momentum_dist1(i,j,k/5,:) = reshape(fp,[1,1,1,Nx]);
            wid(i,j,k) = integr(f.^2.*fp,Nx,deltaf);
            
            
           
            phi = phi_1;
  
        end
        %save(strcat(save_path,'momentum_dist1.mat'),'momentum_dist1')
        save(strcat(save_path,'wid.mat'),'wid')
        
    end
end

% wid = zeros(5,20,25);
% for i=1:5
% for j=1:20
% for k=1:25
% wid(i,j,k) = integr(f.^2.*reshape(momentum_dist1(i,j,k,:),[1,Nx]),Nx,deltaf);
% end
% end
% end
% width = sqrt(reshape(mean(wid,2),[5,25]));