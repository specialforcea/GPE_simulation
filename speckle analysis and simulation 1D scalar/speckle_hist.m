image = spec;
im = image-0;
H = histogram(im);
avg_int = mean(mean(im));
x = H.BinEdges;
y = H.Values;
x = x(1:350)./avg_int;
y = y(1:350);

f = fit(x',log(y)','a*x+b','StartPoint',[-1,10]);
%f = fit(x',y','a*x+b)','StartPoint',[-1,10]);
coe = coeffvalues(f);
plot(f,x,log(y))
xlabel('I/\langleI\rangle')
ylabel('log pixel counts')
text(1,4,strcat('y=',num2str(coe(1)),'x+',num2str(coe(2))))