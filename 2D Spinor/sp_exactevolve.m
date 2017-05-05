function  amp = sp_exactevolve(t,Omega,det)


amp_0 = [1;0;0];



Ham = [0 Omega 0;Omega 0 Omega; 0 Omega det];

[V,D,B] = eig(Ham);

eigens = [D(1,1);D(2,2);D(3,3)];

amp = V*(exp(-1i*eigens*double(t)).*(B'*amp_0));
end




