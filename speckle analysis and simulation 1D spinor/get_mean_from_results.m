folders = [05,1,15,2,25,3,35,4,45,5,55,6,65,7,75,10,11,12];

savepath = 'simulation_results/05162019kick_evolve_with_soc_cluster/';
mp = zeros(18,21,320);
for i=1:18
    result_path = strcat(savepath,num2str(folders(i)));
    load(strcat(result_path,'\mean_p.mat'));
    load(strcat(result_path,'\spin1.mat'));
    mean_p_ = reshape(mean_p(:,:,1,:),[21,20,320]).*spin1 + reshape(mean_p(:,:,2,:),[21,20,320]).*(1-spin1);
    mp(i,:,:) = reshape(mean(mean_p_,2),[1,21,320]);
    gv = (mp(:,:,2:320)-mp(:,:,1:319))/(0.005*Dip_freq);
    gv_ratio = reshape(gv(:,:,319)./gv(:,:,1),[18,21]);
end
