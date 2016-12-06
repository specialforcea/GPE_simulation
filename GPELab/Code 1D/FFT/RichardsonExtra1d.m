%% Computation of a step in time using the full Splitting scheme with Richardson extrapolation
%% INPUTS:
%%          Phi: Initial wave functions in the 1D geometry for the FFT (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          FFTGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FFT context (structure) (see FFTGeometry1D_Var1d.m)
%%          FFTPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FFT context (structure) (see FFTPhysics1D_Var1d.m)
%%          FFTOperators1D: Structure containing the derivative FFT operators (structure) (see FFTOperators1D_Var1d.m)
%% OUTPUTS:
%%          Phi: Wave functions computated with the Richardson extrapolation method on a single step (cell array)

function [Phi , Psi, FFTPsi, flag, relres, iter, resvec] = RichardsonExtra1d(Phi, Psi, FFTPsi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D)

[Rich(1).Phi, Psi_Richtmp, FFTPsi_Richtmp, flag, relres, iter, resvec] = Local_Full_RSP_solution1d(Phi, Psi, FFTPsi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D);

for j = 1:Method.Richardson
    Phi_Richtmp = Phi;
    Psi_Richtmp = Psi;
    FFTPsi_Richtmp = FFTPsi;
    
    Method.Deltat = Method.Deltat/2;
    for k = 1 : 2^j
        [Phi_Richtmp, Psi_Richtmp, FFTPsi_Richtmp] = Local_Full_RSP_solution1d(Phi_Richtmp, Psi_Richtmp, FFTPsi_Richtmp, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D);
    end
    Rich(j+1).Phi = Phi_Richtmp;
end

for j = 1:Method.Richardson
    for k = 1:Method.Richardson + 1  - j
        Rich(k).Phi = Mulc(Addc(Mulc(Rich(k+1).Phi,2^(2*j)),Mulc(Rich(k).Phi,-1)),(2^(2*j)-1)^(-1));
    end
end

Phi = Rich(1).Phi;

for n = 1:Method.Ncomponents
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics1D.Nonlinearity_function_Index{n}
        FFTPhysics1D.Nonlinearity{n,m} = FFTPhysics1D.Nonlinearity_function{n,m}(Phi,FFTGeometry1D.X); % Computing and storing the coupled nonlinearities between components
        Psi{n,m} = 2*FFTPhysics1D.Nonlinearity{n,m} - Psi{n,m}; % Computing the relaxation variable corresponding to the local nonlinearities
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics1D.FFTNonlinearity_function_Index{n}
        FFTPhysics1D.FFTNonlinearity{n,m} = FFTPhysics1D.FFTNonlinearity_function{n,m}(Phi,FFTGeometry1D.X,-1i*FFTOperators1D.Gx); % Computing and storing the coupled non-local nonlinearities between components
        FFTPsi{n,m} = 2*FFTPhysics1D.FFTNonlinearity{n,m} - FFTPsi{n,m}; % Computing the relaxation variable corresponding to the non-local nonlinearities
    end
end