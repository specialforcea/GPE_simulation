%% Computation of a step in time using the Richardson extrapolation
%% INPUTS:
%%          Phi: Initial wave functions (cell array)
%%          Operator: Function with arguments (Phi,Method) to apply the scheme that will be extrapolated (function)
%% OUTPUTS:
%%          Phi: Wave functions computated with the Richardson extrapolation method on a single step (cell array)

function Phi = RichardsonExtra(Phi , Method, Operator)

for j = 1:Method.Richardson
    Phi_Richtmp = Phi;
    Method.Deltat = Method.Deltat/2;
    for k = 1 : 2^j
        Phi_Richtmp = Operator(Phi_Richtmp, Method);
    end
    Phi_Rich{j+1} = Phi_Richtmp;
end

for j = 1:Method.Richardson
    for k = 1:Method.Richardson + 1  - j
        Phi_Rich{k} = Mulc(Addc(Mulc(Phi_Rich{k+1},2^(2+k-1)),Mulc(Phi_Rich{k},-1)),(2^(2+k-1)-1)^(-1));
    end
end

Phi = Phi_Rich{1};