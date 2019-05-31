folders = ["05" "1" "15" "2" "25" "3" "35" "4" "45" "5" "55" "6" "65" "7" "75" "10" "11" "12"];
omegas = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,10,11,12];

savepath = 'simulation_results/05162019kick_evolve_with_soc_cluster/';
mp = zeros(18,21,320);
mm = zeros(18,21,320);
gv = zeros(18,21,320);
for i=1:18
    result_path = strcat(savepath,num2str(folders(i)));
    load(strcat(result_path,'\mean_p.mat'));
    load(strcat(result_path,'\spin1.mat'));
    mean_p_ = reshape(mean_p(:,:,1,:),[21,20,320]).*spin1 + reshape(mean_p(:,:,2,:),[21,20,320]).*(1-spin1);
    mp(i,:,:) = reshape(mean(mean_p_,2),[1,21,320]);
end
gv = (mp(:,:,2:320)-mp(:,:,1:319))/(0.005*Dip_freq);
gv_ratio = reshape(gv(:,:,130)./gv(:,:,1),[18,21]);

% for i=1:18
%     result_path = strcat(savepath,folders(i));
%     load(strcat(result_path,'\mean_m.mat'));
%     load(strcat(result_path,'\spin1.mat'));
%     %mean_m_ = reshape(mean_m(:,:,1,:),[21,20,320]).*spin1 + reshape(mean_m(:,:,2,:),[21,20,320]).*(1-spin1);
%     mean_m_ = reshape(mean_m(:,:,1,:),[21,20,320]);
%     mm(i,:,:) = reshape(mean(mean_m_,2),[1,21,320])./k_R*2*pi;
%     gv(i,:,:) = 2.*mm(i,:,:) + 2 - 8*(mm(i,:,:)+1)./sqrt(16.*(mm(i,:,:)+1).^2 + omegas(i)^2);
% end
% 
% gv_ratio = reshape(gv(:,:,180)./gv(:,:,1),[18,21]);
