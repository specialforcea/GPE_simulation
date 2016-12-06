%% Applying the splitting scheme on a single time step
%% INPUTS:
%%          Phi_in: Initial components' wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          FFTGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FFT context (structure) (see FFTGeometry1D_Var1d.m)
%%          FFTPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FFT context (structure) (see FFTPhysics1D_Var1d.m)
%%          FFTOperators1D: Structure containing the derivative FFT operators (structure) (see FFTOperators1D_Var1d.m)
%% OUTPUT:
%%          Phi: Components' wave functions with the exponential operators applied (cell array)

function [Phi_up] = operator_Splitting1d(Phi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D)
%% Initializing the exponantial operator
MatcPot = cell(Method.Ncomponents); % Initializing the variable that will contain the potential and nonlinear operators

%% Storing starting time
Time = Method.Deltat*(Method.Iterations-1);

if (strcmp(Method.Splitting, 'Lie'))

Time0 = Time;
Time1 = Time + Method.Deltat;
ExpMatcPot = expc(Update_MatPot(Phi)); % Computing the exponential operator for the potential and nonlinear operators

Phi_up = Mulc(Phi,ExpMatcPot);
for n = 1:Method.Ncomponents
    Phi_up{n} = Phi_up{n}/L2_norm1d(Phi_up{n},FFTGeometry1D)*sqrt(Method.NParticles(1)); % Normalization of each wave function
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi_up{n}); 
end

for n = 1:Method.Ncomponents
    FFTPhi{n} = FFTPhi{n}/L2_norm1d(FFTPhi{n},FFTGeometry1D)*sqrt(Method.NParticles(1)); % Normalization of each wave function
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics1D.ExpMatcDiffx,expc(Update_MatDisp)));

for n = 1:Method.Ncomponents
    FFTPhi{n} = FFTPhi{n}/L2_norm1d(FFTPhi{n},FFTGeometry1D)*sqrt(Method.NParticles(1)); % Normalization of each wave function
end
% IFFT
for n = 1:Method.Ncomponents
    Phi_up{n} =ifft(FFTPhi{n}); 
end

for n = 1:Method.Ncomponents
    Phi_up{n} = Phi_up{n}/L2_norm1d(Phi_up{n},FFTGeometry1D)*sqrt(Method.NParticles(1)); % Normalization of each wave function
end

elseif (strcmp(Method.Splitting, 'Strang'))

Time0 = Time;
Time1 = Time + Method.Deltat/2;

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n}); 
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics1D.ExpMatcDiffx,expc(Update_MatDisp)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n}); 
end

Time0 = Time;
Time1 = Time + Method.Deltat;
ExpMatcPot = expc(Update_MatPot(Phi)); % Computing the exponential operator for the potential and nonlinear operators
Phi = Mulc(Phi,ExpMatcPot);

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n}); 
end

Time0 = Time + Method.Deltat/2;
Time1 = Time + Method.Deltat;

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics1D.ExpMatcDiffx,expc(Update_MatDisp)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n}); 
end

elseif (strcmp(Method.Splitting, 'Fourth'))
%% Fourth-order splitting scheme on a single time step
% Computing time-weights
w_1 = 0.33780179798991440851;
w_2 = 0.67560359597982881702;
w_3 = -0.08780179798991440851;
w_4 = -0.85120719795965763405;

% Applying the exponential operators
Time0 = Time;
Time1 = Time + 2*w_1*Method.Deltat;

ExpMatcPot1 = expc(Update_MatPot(Phi)); % Computing the exponential operator for the potential and nonlinear operators
Phi = Mulc(Phi,ExpMatcPot1);

Time0 = Time;
Time1 = Time + 2*w_2*Method.Deltat;

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n}); 
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics1D.ExpMatcDiffx1,expc(Update_MatDisp)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n}); 
end

Time0 = Time + 2*w_1*Method.Deltat;
Time1 = Time + 2*(w_1+w_3)*Method.Deltat;

ExpMatcPot2 = expc(Update_MatPot(Phi)); % Computing the exponential operator for the potential and nonlinear operators
Phi = Mulc(Phi,ExpMatcPot2);

Time0 = Time + 2*w_2*Method.Deltat;
Time1 = Time + 2*(w_2+w_4)*Method.Deltat;

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n}); 
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics1D.ExpMatcDiffx2,expc(Update_MatDisp)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n}); 
end

Time0 = Time + 2*(w_1+w_3)*Method.Deltat;
Time1 = Time + 2*(w_1+2*w_3)*Method.Deltat;

ExpMatcPot2 = expc(Update_MatPot(Phi)); % Computing the exponential operator for the potential and nonlinear operators
Phi = Mulc(Phi,ExpMatcPot2);

Time0 = Time + 2*(w_2+w_4)*Method.Deltat;
Time1 = Time + 2*(2*w_2+w_4)*Method.Deltat;

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n}); 
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics1D.ExpMatcDiffx1,expc(Update_MatDisp)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n}); 
end

Time0 = Time + 2*(w_1+2*w_3)*Method.Deltat;
Time1 = Time + 4*(w_1+w_3)*Method.Deltat;
ExpMatcPot1 = expc(Update_MatPot(Phi)); % Computing the exponential operator for the potential and nonlinear operators
Phi = Mulc(Phi,ExpMatcPot1);

end

function [MatcPot] = Update_MatPot(Phi)

%% Computing the explicit nonlinearity, non-local nonlinearity and the time-depend potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics1D.Nonlinearity_function_Index{n}
        FFTPhysics1D.Nonlinearity{n,m} = FFTPhysics1D.Nonlinearity_function{n,m}(Phi,FFTGeometry1D.X); % Computing and storing the coupled nonlinearities between components
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics1D.FFTNonlinearity_function_Index{n}
        FFTPhysics1D.FFTNonlinearity{n,m} = FFTPhysics1D.FFTNonlinearity_function{n,m}(Phi,FFTGeometry1D.X,-1i*FFTOperators1D.Gx); % Computing and storing the coupled non-local nonlinearities between components
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics1D.TimePotential_function_Index{n}
        % IF the integrated time potential is not defined
        if (isempty(FFTPhysics1D.IntegratedTimePotential_function{n,m}))
            FFTPhysics1D.TimePotential{n,m} = (1/6)*(FFTPhysics1D.TimePotential_function{n,m}(abs(Time0),FFTGeometry1D.X) + 4*FFTPhysics1D.TimePotential_function{n,m}((abs(Time0+Time1))/2,FFTGeometry1D.X) + FFTPhysics1D.TimePotential_function{n,m}(abs(Time1),FFTGeometry1D.X))*(Time1-Time0); % Computing and storing the coupled time-dependent potential between components
        % ELSE if the integrated time potential is defined
        else
            FFTPhysics1D.TimePotential{n,m} = -1i*(FFTPhysics1D.IntegratedTimePotential_function{n,m}(abs(Time1),FFTGeometry1D.X) - FFTPhysics1D.IntegratedTimePotential_function{n,m}(abs(Time0),FFTGeometry1D.X));
        end
    end
    % FOR each component where the stochastic potential is non null
    for m = FFTPhysics1D.StochasticPotential_function_Index{n}
        if (iscell(FFTPhysics1D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics1D.StochasticProcess_function)
                FFTPhysics1D.StochasticProcess{m_noise} = FFTPhysics1D.StochasticProcess_function{m_noise}(abs(Time1),FFTGeometry1D.X)-FFTPhysics1D.StochasticProcess_function{m_noise}(abs(Time0),FFTGeometry1D.X);
            end
        else
            FFTPhysics1D.StochasticProcess = FFTPhysics1D.StochasticProcess_function(abs(Time1),FFTGeometry1D.X)-FFTPhysics1D.StochasticProcess_function(abs(Time0),FFTGeometry1D.X);
        end
        FFTPhysics1D.StochasticPotential{n,m} = -1i*FFTPhysics1D.StochasticPotential_function{n,m}(FFTPhysics1D.StochasticProcess,FFTGeometry1D.X); % Computing and storing the stochastic potential
    end
end

% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        MatcPot{n,m} = 0;
    end
end

%% Constructing the cell matrix using the potential and nonlinear operators
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        MatcPot{n,m} = MatcPot{n,m} + (Time1-Time0)*FFTPhysics1D.Potential{n,m}; % Computing and adding the potential operator
        MatcPot{n,m} = MatcPot{n,m} + (Time1-Time0)*FFTPhysics1D.Beta*FFTPhysics1D.Nonlinearity{n,m}; % Computing and adding the nonlinear operator
        MatcPot{n,m} = MatcPot{n,m} + (Time1-Time0)*FFTPhysics1D.Beta*FFTPhysics1D.FFTNonlinearity{n,m}; % Computing and adding the non-local nonlinear operator
        MatcPot{n,m} = MatcPot{n,m} + FFTPhysics1D.TimePotential{n,m}; % Computing and adding the time-dependent potential operator
        MatcPot{n,m} = MatcPot{n,m} + FFTPhysics1D.StochasticPotential{n,m}; % Computing and adding the stochastic potential operator
    end
end

end

function [MatcDisp] = Update_MatDisp

%% Computing the explicit dispersion and the time-depend dispersion
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the time-dependent dispersion is non null
    for m = FFTPhysics1D.TimeDispersion_function_Index{n}
        % IF the integrated time dispersion is not defined
        if (isempty(FFTPhysics1D.IntegratedTimeDispersion_function{n,m}))
            FFTPhysics1D.TimeDispersion{n,m} = (1/6)*(FFTPhysics1D.TimeDispersion_function{n,m}(abs(Time0),-1i*FFTOperators1D.Gx) + 4*FFTPhysics1D.TimeDispersion_function{n,m}((abs(Time0+Time1))/2,-1i*FFTOperators1D.Gx) + FFTPhysics1D.TimeDispersion_function{n,m}(abs(Time1),-1i*FFTOperators1D.Gx))*(Time1-Time0); % Computing and storing the coupled time-dependent dispersion between components
        % ELSE if the integrated time potential is defined
        else
            FFTPhysics1D.TimeDispersion{n,m} = -1i*(FFTPhysics1D.IntegratedTimeDispersion_function{n,m}(abs(Time1),-1i*FFTOperators1D.Gx) - FFTPhysics1D.IntegratedTimeDispersion_function{n,m}(abs(Time0),-1i*FFTOperators1D.Gx));
        end
    end
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics1D.StochasticDispersion_function_Index{n}
        if (iscell(FFTPhysics1D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics1D.StochasticProcess_function)
                FFTPhysics1D.StochasticProcess{m_noise} = FFTPhysics1D.StochasticProcess_function{m_noise}(abs(Time1),FFTGeometry1D.X)-FFTPhysics1D.StochasticProcess_function{m_noise}(abs(Time0),FFTGeometry1D.X);
            end
        else
            FFTPhysics1D.StochasticProcess = FFTPhysics1D.StochasticProcess_function(abs(Time1),FFTGeometry1D.X)-FFTPhysics1D.StochasticProcess_function(abs(Time0),FFTGeometry1D.X);
        end
        FFTPhysics1D.StochasticDispersion{n,m} = -1i*FFTPhysics1D.StochasticDispersion_function{n,m}(FFTPhysics1D.StochasticProcess,-1i*FFTOperators1D.Gx); % Computing and storing the stochastic potential
    end
end

% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        MatcDisp{n,m} = 0;
    end
end

%% Constructing the cell matrix using the potential and nonlinear operators
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        MatcDisp{n,m} = MatcDisp{n,m} + FFTPhysics1D.StochasticDispersion{n,m} + FFTPhysics1D.TimeDispersion{n,m}; % Computing and adding the potential operator
    end
end

end
end