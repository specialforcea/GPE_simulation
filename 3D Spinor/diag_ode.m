mat = [9*A+2*B B 0 0 0 0 0;B 4*A+2*B B 0 0 0 0;0 B A+2*B B 0 0 0;0 0 B 2*B B 0 0;0 0 0 B A+2*B B 0;0 0 0 0 B 4*A+2*B B;0 0 0 0 0 B 9*A+2*B];
                                                                                                                             
s = size(t_ode,1);
c = zeros(s,7);
for j = 1:s
   [c(j,1),c(j,2),c(j,3),c(j,4),c(j,5),c(j,6),c(j,7)] = exact_ode(mat,t_ode(j));
end

