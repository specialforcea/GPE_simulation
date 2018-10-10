% for i=1:4
% for j=1:6
% fpath = strcat('evolve under speckle/',num2str(i),num2str(j),'.png');
% plot(X.*a_0./1e-6,reshape(phi_all(i,j,:),[1,Nx]));
% 
% tit = strcat('cor len 1/',num2str(2^i),' intensity ', num2str(intensity(j)));
% title(tit);
% xlabel('x/um');
% ylabel('|phi|^2');
% saveas(gcf,fpath);
% end
% end
x = X.*deltax.*a_0./1e-6;
% for i=1:4
%     phi1 = reshape(phi_all(i,1,:),[1,Nx]);
%    
%     phi3 = reshape(phi_all(i,3,:),[1,Nx]);
%     phi4 = reshape(phi_all(i,4,:),[1,Nx]);
%     
%     phi6 = reshape(phi_all(i,6,:),[1,Nx]);
%     plot(x,phi1,x,phi3,x,phi4,x,phi6);
%     
%     tit = strcat('cor len ',num2str(2^(i-1)),' um different intensity');
%     title(tit);
%     xlabel('x/um');
%     ylabel('|phi|^2');
%     legend('int 1/4','int 1/2','int 1','int 3');
%     fpath = strcat('evolve under speckle/','cor len ',num2str(2^(i-1)),' um','.png');
%     saveas(gcf,fpath);
%     
% end

for j=[1,3,4,6]
    phi1 = reshape(phi_all(1,j,:),[1,Nx]);
   
    phi3 = reshape(phi_all(2,j,:),[1,Nx]);
    phi4 = reshape(phi_all(3,j,:),[1,Nx]);
    
    phi6 = reshape(phi_all(4,j,:),[1,Nx]);
    plot(x,phi1,x,phi3,x,phi4,x,phi6);
    
    tit = strcat('intensity ',num2str(intensity(j)),' different cor len');
    title(tit);
    xlabel('x/um');
    ylabel('|phi|^2');
    legend('cor len 1 um','cor len 2 um','cor len 4 um','cor len 8 um');
    fpath = strcat('evolve under speckle/','intensity ',num2str(intensity(j)),'.png');
    saveas(gcf,fpath);
    
end