function [Potential] = Gaussian_potential(X,Y,w,W,d)

Potential = W*exp(-4*((X-d(1)).^2 + (Y-d(2)).^2)/w^2);