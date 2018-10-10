xcor_len = zeros(1,4);
ycor_len = zeros(1,4);
for i=1:4
    filepath = strcat('speckle bench test data/9/',num2str(i),'.png');
    image = imread(filepath);
    image = double(image);
    
    [xcor,ycor] = correlation_length_of_2d_image(image);
    [~,nx] = size(xcor);
    [~,ny] = size(ycor);
    
    x = (1:1:100).*4.8;
    y = (1:1:100).*4.8;
    
    l = find(xcor<=0.135,1)*4.8;
    if isempty(l)
        xcor_len(1,i) = 0;
    else
        xcor_len(1,i) = l;
    end
    
    l = find(ycor<=0.135,1)*4.8;
    if isempty(l)
        ycor_len(1,i) = 0;
    else
        ycor_len(1,i) = l;
    end
 
    
    plot(x,xcor(1:100));
    ti = strcat('x correlation for image ',num2str(i),'.png');
    title(ti);
    xlabel('delta_x/um');
    ylabel('correlation')
    saveas(gcf,strcat('correlation images/9/',ti))
    
    
    
    plot(y,ycor(1:100));
    ti = strcat('y correlation for image ',num2str(i),'.png');
    title(ti);
    xlabel('delta_y/um');
    ylabel('correlation')
    saveas(gcf,strcat('correlation images/9/',ti))
    
end