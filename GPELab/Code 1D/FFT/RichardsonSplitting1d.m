%% Computation of a step in time using the full Splitting scheme with Richardson extrapolation
%% INPUTS:
%%          Phi: Initial wave functions in the 1D geometry for the FFT (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          FFTGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FFT context (structure) (see FFTGeometry1D_Var1d.m)
%%          FFTPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FFT context (structure) (see FFTPhysics1D_Var1d.m)
%%          FFTOperators1D: Structure containing the derivative FFT operators (structure) (see FFTOperators1D_Var1d.m)
%% OUTPUTS:
%%          Phi: Wave functions computated with the Richardson extrapolation method on a single step (cell array)

function Phi = RichardsonSplitting1d(Phi,Method,FFTGeometry1D,FFTPhysics1D,FFTOperators1D)

Rich(1).Phi = Local_Splitting_solution1d(Phi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D);

for j = 1:Method.Richardson
    Phi_Richtmp = Phi;
    Method.Deltat = Method.Deltat/2;
    for k = 1 : 2^j
        Phi_Richtmp = Local_Splitting_solution1d(Phi_Richtmp, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D);
    end
    Rich(j+1).Phi = Phi_Richtmp;
end

if (strcmp(Method.Splitting,'Lie'))
    l = 1;
elseif (strcmp(Method.Splitting,'Strang'))
    l = 2;
elseif (strcmp(Method.Splitting,'Fourth'))
    l = 4;
end

for j = 1:Method.Richardson
    for k = 1:(Method.Richardson - j + 1)
        Rich(k).Phi = Mulc(Addc(Mulc(Rich(k+1).Phi,2^(l*j)),Mulc(Rich(k).Phi,-1)),(2^(l*j)-1)^(-1));
    end
end

Phi = Rich(1).Phi;