ho = [0.0841;0;0];
shit = zeros(1,40000);
for pp = 1:40000
    ho = MatV(:,:,25444)*(exp(-eigens(:,25444).*Deltat*1i).*(MatB(:,:,25444)*ho));
    shit(pp) = ho(1);
end
x = 1:40000;
