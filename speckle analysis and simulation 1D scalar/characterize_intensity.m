image = spec;
im = image-0;
[Ny, Nx] = size(image);
%im = double(im-min(min(im)));

V_average = mean(mean(im));
max_im = floor((max(max(im))));
x = (0:1:max_im)/V_average;
y = zeros(1,max_im+1);
for i=Ny
for j=Nx
y(1,floor(im(i,j)) + 1) = y(1,floor(im(i,j))+1) + 1;
end
end
% y = y/sum(y);
% max_ind = y==max(y);
% y = y*exp(-x(max_ind))/max(y);
% y = log(y);
% z = -x;
plot(x,log(y));
xlabel('I/<I>');
ylabel('log[P(I)]');
legend('data','exp(-x)');
title('speckle intensity distibution');
