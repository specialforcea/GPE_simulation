function V_average = characterize_intensity(image,dark)
im = image-dark;
im = double(im-min(min(im)));

V_average = mean(mean(im));
max_im = max(max(im));
x = (0:1:max_im)/V_average/0.30;
y = zeros(1,max_im+1);
for i=1:1024
for j=1:1280
y(im(i,j) + 1) = y(im(i,j)+1) + 1;
end
end
y = y/sum(y);
max_ind = y==max(y);
y = y*exp(-x(max_ind))/max(y);
y = log(y);
z = -x;
plot(x,y,x,z);
xlabel('I/<I>');
ylabel('log[P(I)]');
legend('data','exp(-x)');
title('speckle intensity distibution');
