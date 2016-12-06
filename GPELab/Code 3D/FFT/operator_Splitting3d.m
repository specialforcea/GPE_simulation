%% Computing the splitting scheme on a single step time
%% INPUTS:
%%          Phi: Initial components' wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          FFTGeometry3D: Structure containing variables concerning the geometry of the problem in 3D in the FFT context (structure) (see FFTGeometry3D_Var3d.m)
%%          FFTPhysics3D: Structure containing variables concerning the physics of the problem in 3D in the FFT context (structure) (see FFTPhysics3D_Var3d.m)
%%          FFTOperators3D: Structure containing the derivative FFT operators (structure) (see FFTOperators3D_Var3d.m)
%% OUTPUT:
%%          Phi: Components' wave functions with the operators applied (cell array)

function [Phi] = operator_Splitting3d(Phi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D)
%% Initializing the exponantial operator
MatcPot = cell(Method.Ncomponents); % Initializing the variable that will contain the potential and nonlinear operators

%% Storing starting time
Time = Method.Deltat*(Method.Iterations-1);

if (strcmp(Method.Splitting, 'Lie'))
    
Time0 = Time;
Time1 = Time + Method.Deltat;
ExpMatcPot = expc(Update_MatPot(Phi)); % Computing the exponential operator for the potential and nonlinear operators
Phi = Mulc(Phi,ExpMatcPot);

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],2); 
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffx,expc(Update_MatDispx)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],2);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],3); 
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffz,expc(Update_MatDispz)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],3);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],1); 
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffy,expc(Update_MatDispy)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],1);  
end

elseif (strcmp(Method.Splitting, 'Strang'))

Time0 = Time;
Time1 = Time + Method.Deltat/2;

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],2);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffx,expc(Update_MatDispx)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],2);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],3);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffz,expc(Update_MatDispz)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],3);  
end


% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],1);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffy,expc(Update_MatDispy)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],1);  
end

Time0 = Time;
Time1 = Time + Method.Deltat;
ExpMatcPot = expc(Update_MatPot(Phi)); % Computing the exponential operator for the potential and nonlinear operators
Phi = Mulc(Phi,ExpMatcPot);
Time0 = Time + Method.Deltat/2;
Time1 = Time + Method.Deltat;

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],1);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffy,expc(Update_MatDispy)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],1);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],3);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffz,expc(Update_MatDispz)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],3);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],2);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffx,expc(Update_MatDispx)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],2);  
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

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],2);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffx1,expc(Update_MatDispx)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],2);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],3);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffz1,expc(Update_MatDispz)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],3);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],1);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffy1,expc(Update_MatDispy)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],1);  
end


Time0 = Time;
Time1 = Time + 2*w_2*Method.Deltat;
ExpMatcPot = expc(Update_MatPot(Phi)); % Computing the exponential operator for the potential and nonlinear operators
Phi = Mulc(Phi,ExpMatcPot);

Time0 = Time + 2*w_1*Method.Deltat;
Time1 = Time + 2*(w_1+w_3)*Method.Deltat;

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],1);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffy2,expc(Update_MatDispy)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],1);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],3);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffz2,expc(Update_MatDispz)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],3);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],2);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffx2,expc(Update_MatDispx)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],2);  
end


Time0 = Time + 2*w_2*Method.Deltat;
Time1 = Time + 2*(w_2+w_4)*Method.Deltat;
ExpMatcPot = expc(Update_MatPot(Phi)); % Computing the exponential operator for the potential and nonlinear operators
Phi = Mulc(Phi,ExpMatcPot);

Time0 = Time + 2*(w_1+w_3)*Method.Deltat;
Time1 = Time + 2*(2*w_1+w_3)*Method.Deltat;

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],2);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffx2,expc(Update_MatDispx)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],2);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],3);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffz2,expc(Update_MatDispz)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],3);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],1);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffy2,expc(Update_MatDispy)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],1);  
end

Time0 = Time + 2*(w_2+w_4)*Method.Deltat;
Time1 = Time + 2*(2*w_2+w_4)*Method.Deltat;
ExpMatcPot = expc(Update_MatPot(Phi)); % Computing the exponential operator for the potential and nonlinear operators
Phi = Mulc(Phi,ExpMatcPot);

Time0 = Time + 2*(2*w_1+w_3)*Method.Deltat;
Time1 = Time + 2*(2*w_1+2*w_3)*Method.Deltat;

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],1);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffy1,expc(Update_MatDispy)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],1);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],3);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffz1,expc(Update_MatDispz)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],3);  
end


% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],2);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics3D.ExpMatcDiffx1,expc(Update_MatDispx)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],2);  
end

end

function [MatcPot] = Update_MatPot(Phi)
%% Computing the explicit nonlinearity, non-local nonlinearity and the time-depend potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics3D.Nonlinearity_function_Index{n}
        FFTPhysics3D.Nonlinearity{n,m} = FFTPhysics3D.Nonlinearity_function{n,m}(Phi,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing and storing the coupled nonlinearities between components
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics3D.FFTNonlinearity_function_Index{n}
        FFTPhysics3D.FFTNonlinearity{n,m} = FFTPhysics3D.FFTNonlinearity_function{n,m}(Phi,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z,-1i*FFTOperators3D.Gx,-1i*FFTOperators3D.Gy,-1i*FFTOperators3D.Gz); % Computing and storing the coupled non-local nonlinearities between components
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics3D.TimePotential_function_Index{n}
        % IF the integrated time potential is not defined
        if (isempty(FFTPhysics3D.IntegratedTimePotential_function{n,m}))
            FFTPhysics3D.TimePotential{n,m} = ((Time1-Time0)/6)*(FFTPhysics3D.TimePotential_function{n,m}(abs(Time1),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z) + FFTPhysics3D.TimePotential_function{n,m}(abs(Time0),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z) +4*FFTPhysics3D.TimePotential_function{n,m}((abs(Time1+Time0))/2,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)); % Computing and storing the coupled time-dependent potential between components
        % ELSE if the integrated time potential is defined
        else
            FFTPhysics3D.TimePotential{n,m} = -1i*(FFTPhysics3D.IntegratedTimePotential_function{n,m}(abs(Time1),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z) - FFTPhysics3D.IntegratedTimePotential_function{n,m}(abs(Time0),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z));
        end
        
    end
    % FOR each component where the stochastic potential is non null
    for m = FFTPhysics3D.StochasticPotential_function_Index{n}
        if (iscell(FFTPhysics3D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics3D.StochasticProcess_function)
                FFTPhysics3D.StochasticProcess{m_noise} = FFTPhysics3D.StochasticProcess_function{m_noise}(abs(Time1),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)-FFTPhysics3D.StochasticProcess_function{m_noise}(abs(Time0),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z);
            end
        else
            FFTPhysics3D.StochasticProcess = FFTPhysics3D.StochasticProcess_function(abs(Time1),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)-FFTPhysics3D.StochasticProcess_function(abs(Time0),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z);
        end
        FFTPhysics3D.StochasticPotential{n,m} = -1i*FFTPhysics3D.StochasticPotential_function{n,m}(FFTPhysics3D.StochasticProcess,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing and storing the stochastic potential
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
        MatcPot{n,m} = MatcPot{n,m} + (Time1-Time0)*FFTPhysics3D.Potential{n,m}; % Computing and adding the potential operator
        MatcPot{n,m} = MatcPot{n,m} + (Time1-Time0)*FFTPhysics3D.Beta*FFTPhysics3D.Nonlinearity{n,m}; % Computing and adding the nonlinear operator
        MatcPot{n,m} = MatcPot{n,m} + (Time1-Time0)*FFTPhysics3D.Beta*FFTPhysics3D.FFTNonlinearity{n,m}; % Computing and adding the non-local nonlinear operator
        MatcPot{n,m} = MatcPot{n,m} + FFTPhysics3D.TimePotential{n,m}; % Computing and adding the time-dependent potential operator
        MatcPot{n,m} = MatcPot{n,m} + FFTPhysics3D.StochasticPotential{n,m}; % Computing and adding the stochastic potential operator
    end
end
end
function [MatcDisp] = Update_MatDispx

%% Computing the explicit nonlinearity, non-local nonlinearity and the time-depend potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the time-dependent dispersion is non null
    for m = FFTPhysics3D.TimeDispersion_function_Index{n}
        % IF the integrated time dispersion is not defined
        if (isempty(FFTPhysics3D.IntegratedTimeDispersion_function{n,m}))
            FFTPhysics3D.TimeDispersion{n,m} = ((Time1-Time0)/6)*(FFTPhysics3D.TimeDispersion_function{n,m}(abs(Time0),-1i*FFTOperators3D.Gx,0,0) + 4*FFTPhysics3D.TimeDispersion_function{n,m}((abs(Time0 + Time1))/2,-1i*FFTOperators3D.Gx,0,0) + FFTPhysics3D.TimeDispersion_function{n,m}(abs(Time1),-1i*FFTOperators3D.Gx,0,0)); % Computing and storing the coupled time-dependent dispersion between components
        % ELSE if the integrated time dispersion is defined
        else
            FFTPhysics3D.TimeDispersion{n,m} = -1i*(FFTPhysics3D.IntegratedTimeDispersion_function{n,m}(abs(Time1),-1i*FFTOperators3D.Gx,0,0) - FFTPhysics3D.IntegratedTimeDispersion_function{n,m}(abs(Time0),-1i*FFTOperators3D.Gx,0,0));
        end
    end
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics3D.StochasticDispersion_function_Index{n}
        if (iscell(FFTPhysics3D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics3D.StochasticProcess_function)
                FFTPhysics3D.StochasticProcess{m_noise} = FFTPhysics3D.StochasticProcess_function{m_noise}(abs(Time1),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)-FFTPhysics3D.StochasticProcess_function{m_noise}(abs(Time0),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z);
            end
        else
            FFTPhysics3D.StochasticProcess = FFTPhysics3D.StochasticProcess_function(abs(Time1),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)-FFTPhysics3D.StochasticProcess_function(abs(Time0),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z);
        end
        FFTPhysics3D.StochasticDispersion{n,m} = -1i*FFTPhysics3D.StochasticDispersion_function{n,m}(FFTPhysics3D.StochasticProcess,-1i*FFTOperators3D.Gx,0,0); % Computing and storing the stochastic potential
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
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics3D.StochasticDispersion_function_Index{n}
        MatcDisp{n,m} = MatcDisp{n,m} + FFTPhysics3D.StochasticDispersion{n,m}; % Computing and adding the stochastic dispersion operator
    end
    % FOR each component where the time-dependent dispersion is non null
    for m = FFTPhysics3D.TimeDispersion_function_Index{n}
        MatcDisp{n,m} = MatcDisp{n,m} + FFTPhysics3D.TimeDispersion{n,m}; % Computing and adding the time-dependent dispersion operator
    end
end

end
function [MatcDisp] = Update_MatDispy

%% Computing the explicit nonlinearity, non-local nonlinearity and the time-depend potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the time-dependent dispersion is non null
    for m = FFTPhysics3D.TimeDispersion_function_Index{n}
        % IF the integrated time dispersion is not defined
        if (isempty(FFTPhysics3D.IntegratedTimeDispersion_function{n,m}))
            FFTPhysics3D.TimeDispersion{n,m} = ((Time1-Time0)/6)*(FFTPhysics3D.TimeDispersion_function{n,m}(abs(Time0),0,-1i*FFTOperators3D.Gy,0) + 4*FFTPhysics3D.TimeDispersion_function{n,m}((abs(Time0 + Time1))/2,0,-1i*FFTOperators3D.Gy,0) + FFTPhysics3D.TimeDispersion_function{n,m}(abs(Time1),0,-1i*FFTOperators3D.Gy,0)); % Computing and storing the coupled time-dependent dispersion between components
        % ELSE if the integrated time dispersion is defined
        else
            FFTPhysics3D.TimeDispersion{n,m} = -1i*(FFTPhysics3D.IntegratedTimeDispersion_function{n,m}(abs(Time1),0,-1i*FFTOperators3D.Gy,0) - FFTPhysics3D.IntegratedTimeDispersion_function{n,m}(abs(Time0),0,-1i*FFTOperators3D.Gy,0));
        end
    end
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics3D.StochasticDispersion_function_Index{n}
        if (iscell(FFTPhysics3D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics3D.StochasticProcess_function)
                FFTPhysics3D.StochasticProcess{m_noise} = FFTPhysics3D.StochasticProcess_function{m_noise}(abs(Time1),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)-FFTPhysics3D.StochasticProcess_function{m_noise}(abs(Time0),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z);
            end
        else
            FFTPhysics3D.StochasticProcess = FFTPhysics3D.StochasticProcess_function(abs(Time1),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)-FFTPhysics3D.StochasticProcess_function(abs(Time0),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z);
        end
        FFTPhysics3D.StochasticDispersion{n,m} = -1i*FFTPhysics3D.StochasticDispersion_function{n,m}(FFTPhysics3D.StochasticProcess,0,-1i*FFTOperators3D.Gy,0); % Computing and storing the stochastic potential
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
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics3D.StochasticDispersion_function_Index{n}
        MatcDisp{n,m} = MatcDisp{n,m} + FFTPhysics3D.StochasticDispersion{n,m}; % Computing and adding the stochastic dispersion operator
    end
        % FOR each component where the time-dependent dispersion is non null
    for m = FFTPhysics3D.TimeDispersion_function_Index{n}
        MatcDisp{n,m} = MatcDisp{n,m} + FFTPhysics3D.TimeDispersion{n,m}; % Computing and adding the time-dependent dispersion operator
    end
end

end

function [MatcDisp] = Update_MatDispz

%% Computing the explicit nonlinearity, non-local nonlinearity and the time-depend potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the time-dependent dispersion is non null
    for m = FFTPhysics3D.TimeDispersion_function_Index{n}
        % IF the integrated time dispersion is not defined
        if (isempty(FFTPhysics3D.IntegratedTimeDispersion_function{n,m}))
            FFTPhysics3D.TimeDispersion{n,m} = ((Time1-Time0)/6)*(FFTPhysics3D.TimeDispersion_function{n,m}(abs(Time0),0,0,-1i*FFTOperators3D.Gz) + 4*FFTPhysics3D.TimeDispersion_function{n,m}((abs(Time0 + Time1))/2,0,0,-1i*FFTOperators3D.Gz) + FFTPhysics3D.TimeDispersion_function{n,m}(abs(Time1),0,0,-1i*FFTOperators3D.Gz)); % Computing and storing the coupled time-dependent dispersion between components
        % ELSE if the integrated time dispersion is defined
        else
            FFTPhysics3D.TimeDispersion{n,m} = -1i*(FFTPhysics3D.IntegratedTimeDispersion_function{n,m}(abs(Time1),0,0,-1i*FFTOperators3D.Gz) - FFTPhysics3D.IntegratedTimeDispersion_function{n,m}(abs(Time0),0,0,-1i*FFTOperators3D.Gz));
        end
    end
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics3D.StochasticDispersion_function_Index{n}
        if (iscell(FFTPhysics3D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics3D.StochasticProcess_function)
                FFTPhysics3D.StochasticProcess{m_noise} = FFTPhysics3D.StochasticProcess_function{m_noise}(abs(Time1),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)-FFTPhysics3D.StochasticProcess_function{m_noise}(abs(Time0),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z);
            end
        else
            FFTPhysics3D.StochasticProcess = FFTPhysics3D.StochasticProcess_function(abs(Time1),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)-FFTPhysics3D.StochasticProcess_function(abs(Time0),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z);
        end
        FFTPhysics3D.StochasticDispersion{n,m} = -1i*FFTPhysics3D.StochasticDispersion_function{n,m}(FFTPhysics3D.StochasticProcess,0,0,-1i*FFTOperators3D.Gz); % Computing and storing the stochastic potential
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
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics3D.StochasticDispersion_function_Index{n}
        MatcDisp{n,m} = MatcDisp{n,m} + FFTPhysics3D.StochasticDispersion{n,m}; % Computing and adding the stochastic dispersion operator
    end
    % FOR each component where the time-dependent dispersion is non null
    for m = FFTPhysics3D.TimeDispersion_function_Index{n}
        MatcDisp{n,m} = MatcDisp{n,m} + FFTPhysics3D.TimeDispersion{n,m}; % Computing and adding the time-dependent dispersion operator
    end
end

end
end