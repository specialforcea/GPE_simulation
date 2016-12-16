mat = [4*A+2*B B 0 0 0;B A+2*B B 0 0;0 B 2*B B 0;0 0 B A+2*B B;0 0 0 B 4*A+2*B];
[right,eigen,left] = eig(mat);
s = size(t_ode,1);
c = zeros(s,5);
for j = 1:s
   [c(j,1),c(j,2),c(j,3),c(j,4),c(j,5)] = exact_ode(mat,t_ode(j));
end

