rho = @(x)(( acos(min(1,x./3))- x./3.*sqrt(1-(min(1,x./3)).^2))*2/pi);
N = 1001;
k = linspace(-2,2,N);
T = zeros(N,N);
E_window = 0.8;
initial = zeros(N,1);
initial(130) = 1;
for i = 1:N
    for j=1:N
        if abs(k(j)^2-k(i)^2)>E_window
            T(i,j) = 0;
        else
            T(i,j) = rho(abs(k(j)-k(i)));
        end
    end
end

% norm = repmat(sum(T,1),N,1);
% T = T./norm;
% [V,D] = eig(T);
diff = 1;
state = initial;
while diff>0.0001
    state1 = T*state;
    state1 = state1./sum(state1);
    diff = max(abs(state1-state));
    state = state1;
end
