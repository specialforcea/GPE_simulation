files = dir('speckle bench test data/for hist/*.png');
intercept = zeros(1,20);
j = 1;
for file=files'
    speckle = imread(strcat(file.folder,'/',file.name));
    spec = speckle - 20;
    image = spec;
    im = image-0;
    H = histogram(im);
    avg_int = mean(mean(im));
    x = H.BinEdges;
    y = H.Values;
    x = (x(2:end)+x(1:end-1))./2./avg_int;

    x = x(25:end);
    y = y(25:end);
    for i=1:size(y,2)
        if y(i)==0
            x(i) = 0;
        end
    end

    x = x(x~=0);
    y = y(y~=0);

    x = x(1:end-10);
    y = y(1:end-10);

    f = fit(x',log(y)','a*x+b','StartPoint',[-1,10]);
    %f = fit(x',y','a*x+b)','StartPoint',[-1,10]);
    coe = coeffvalues(f);
    plot(f,x,log(y))
    xlabel('I/\langleI\rangle')
    ylabel('log pixel counts')
    text(1,4,strcat('y=',num2str(coe(1)),'x+',num2str(coe(2))))
    intercept(j) = coe(1);
    j = j+1;
end