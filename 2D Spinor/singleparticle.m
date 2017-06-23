evo_time = 0.2;
Deltat_count = 1e-4;
Deltat = 1e-5;
evoN = int16(evo_time/Deltat_count);
real_Omega = [1 2 3 4 5 6]*2.34*1e4;
det = [2 1 0 -1 -2];
fprsp = zeros(evoN,3,6,5);
Omega = real_Omega/Dip_freq;
for yy = 1:6
    for yyy = 1:5
        t = 0;
        for i = 1:evoN
            
            fprsp(i,:,yy,yyy) = sq(sp_exactevolve(t,Omega(yy),det(yyy)*detuning));
            t = t + Deltat_count;
        end
    end
end