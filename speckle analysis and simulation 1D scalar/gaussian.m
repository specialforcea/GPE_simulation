function y = gaussian(x,u,sig)
y = 1/sqrt(2*pi*sig^2)* exp(-(x-u).^2/(2*sig^2));
