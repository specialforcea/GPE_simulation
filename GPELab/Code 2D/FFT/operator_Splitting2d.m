%% Applying the splitting scheme on a single time step
%% INPUTS:
%%          Phi_in: Initial components' wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          FFTGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%%          FFTPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FFT context (structure) (see FFTPhysics2D_Var2d.m)
%%          FFTOperators2D: Structure containing the derivative FFT operators (structure) (see FFTOperators2D_Var2d.m)
%% OUTPUT:
%%          Phi: Components' wave functions with the exponential operators applied (cell array)

function [Phi_up] = operator_Splitting2d(Phi, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D)

%% Storing starting time
Time = Method.Deltat*(Method.Iterations-1);

if (strcmp(Method.Splitting, 'Lie'))
%% Lie splitting scheme on a single time step

Time0 = Time;
Time1 = Time + Method.Deltat;
MatcPot = Update_MatPot(Phi);
ExpMatcPot = expc(MatcPot); % Computing the exponential operator for the potential and nonlinear operators

Phi_up = Mulc(Phi,ExpMatcPot);
for n = 1:Method.Ncomponents
   Phi_up{n} = Phi_up{n}/L2_norm2d(Phi_up{n},FFTGeometry2D)*sqrt(Method.NParticles(n)); % Normalization of each wave function
        
end


% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi_up{n},[],2); 
end

for n = 1:Method.Ncomponents
   FFTPhi{n} = FFTPhi{n}/L2_norm2d(FFTPhi{n},FFTGeometry2D)*sqrt(Method.NParticles(n)); % Normalization of each wave function
        
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffx,expc(Update_MatDispx)));
for n = 1:Method.Ncomponents
   FFTPhi{n} = FFTPhi{n}/L2_norm2d(FFTPhi{n},FFTGeometry2D)*sqrt(Method.NParticles(n)); % Normalization of each wave function
        
end


% IFFT
for n = 1:Method.Ncomponents
    Phi_up{n} =ifft(FFTPhi{n},[],2);  
end



for n = 1:Method.Ncomponents
   Phi_up{n} = Phi_up{n}/L2_norm2d(Phi_up{n},FFTGeometry2D)*sqrt(Method.NParticles(n)); % Normalization of each wave function
        
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi_up{n},[],1); 
end

for n = 1:Method.Ncomponents
   FFTPhi{n} = FFTPhi{n}/L2_norm2d(FFTPhi{n},FFTGeometry2D)*sqrt(Method.NParticles(n)); % Normalization of each wave function
        
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffy,expc(Update_MatDispy)));
for n = 1:Method.Ncomponents
   FFTPhi{n} = FFTPhi{n}/L2_norm2d(FFTPhi{n},FFTGeometry2D)*sqrt(Method.NParticles(n)); % Normalization of each wave function
        
end


% IFFT
for n = 1:Method.Ncomponents
    Phi_up{n} =ifft(FFTPhi{n},[],1);  
end

for n = 1:Method.Ncomponents
   Phi_up{n} = Phi_up{n}/L2_norm2d(Phi_up{n},FFTGeometry2D)*sqrt(Method.NParticles(n)); % Normalization of each wave function
        
end



elseif (strcmp(Method.Splitting, 'Strang'))

Time0 = Time;
Time1 = Time + Method.Deltat/2;    

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],2);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffx,expc(Update_MatDispx)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],2);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],1);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffy,expc(Update_MatDispy)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],1);  
end


Time0 = Time;
Time1 = Time + Method.Deltat;
MatcPot = Update_MatPot(Phi);
ExpMatcPot = expc(MatcPot); % Computing the exponential operator for the potential and nonlinear operators


Phi = Mulc(Phi,ExpMatcPot);


Time0 = Time + Method.Deltat/2;
Time1 = Time + Method.Deltat;

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],1);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffy,expc(Update_MatDispy)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],1);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],2);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffx,expc(Update_MatDispx)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],2);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft2(Phi{n}); 
end
    
elseif (strcmp(Method.Splitting, 'Fourth'))

% Computing time-weights
w_1 = 0.33780179798991440851;
w_3 = -0.08780179798991440851;
w_2 = 0.67560359597982881702;
w_4 = -0.85120719795965763405;

Time0 = Time;
Time1 = Time + 2+w_1*Method.Deltat;

% Applying the exponential operators

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],2);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffx1,expc(Update_MatDispx)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],2);
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],1);
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffy1,expc(Update_MatDispy)));

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
    FFTPhi{n} =fft(Phi{n},[],2);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffx2,expc(Update_MatDispx)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],2);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],1);
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffy2,expc(Update_MatDispy)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],1);  
end

Time0 = Time + 2*w_2*Method.Deltat;
Time1 = Time + 2*(w_2+w_4)*Method.Deltat;
ExpMatcPot = expc(Update_MatPot(Phi)); % Computing the exponential operator for the potential and nonlinear operators
Phi = Mulc(Phi,ExpMatcPot);

Time0 = Time + 2*(w_1+w_3)*Method.Deltat;
Time1 = Time + 2*(2*w_1+w_3)*Method.Deltat;

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],1);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffy2,expc(Update_MatDispy)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],1);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],2);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffx2,expc(Update_MatDispx)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],2);  
end

Time0 = Time + 2*(w_2+w_4)*Method.Deltat;
Time1 = Time + 2*(2*w_2+w_4)*Method.Deltat;
ExpMatcPot = expc(Update_MatPot(Phi)); % Computing the exponential operator for the potential and nonlinear operators
Phi = Mulc(Phi,ExpMatcPot);

Time0 = Time + 2*(2*w_1+w_3)*Method.Deltat;
Time1 = Time + 4*(w_1+w_3)*Method.Deltat;

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],1);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffy1,expc(Update_MatDispy)));

% IFFT
for n = 1:Method.Ncomponents
    Phi{n} =ifft(FFTPhi{n},[],1);  
end

% FFT
for n = 1:Method.Ncomponents
    FFTPhi{n} =fft(Phi{n},[],2);  
end

FFTPhi = Mulc(FFTPhi,Mulc(FFTPhysics2D.ExpMatcDiffx1,expc(Update_MatDispx)));

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
    for m = FFTPhysics2D.Nonlinearity_function_Index{n}
        FFTPhysics2D.Nonlinearity{n,m} = FFTPhysics2D.Nonlinearity_function{n,m}(Phi,FFTGeometry2D.X,FFTGeometry2D.Y); % Computing and storing the coupled nonlinearities between components
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics2D.FFTNonlinearity_function_Index{n}
        FFTPhysics2D.FFTNonlinearity{n,m} = FFTPhysics2D.FFTNonlinearity_function{n,m}(Phi,FFTGeometry2D.X,FFTGeometry2D.Y,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Computing and storing the coupled non-local nonlinearities between components
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics2D.TimePotential_function_Index{n}
        % IF the integrated time potential is not defined
        if (isempty(FFTPhysics2D.IntegratedTimePotential_function{n,m}))
            FFTPhysics2D.TimePotential{n,m} =((Time1-Time0)/6)*(FFTPhysics2D.TimePotential_function{n,m}(abs(Time0),FFTGeometry2D.X,FFTGeometry2D.Y) + 4*FFTPhysics2D.TimePotential_function{n,m}((abs(Time0 + Time1))/2,FFTGeometry2D.X,FFTGeometry2D.Y) + FFTPhysics2D.TimePotential_function{n,m}(abs(Time1),FFTGeometry2D.X,FFTGeometry2D.Y)); % Computing and storing the coupled time-dependent potential between components
        % ELSE if the integrated time potential is defined
        else
            FFTPhysics2D.TimePotential{n,m} = -1i*(FFTPhysics2D.IntegratedTimePotential_function{n,m}(abs(Time1),FFTGeometry2D.X,FFTGeometry2D.Y) - FFTPhysics2D.IntegratedTimePotential_function{n,m}(abs(Time0),FFTGeometry2D.X,FFTGeometry2D.Y));
        end
    end
    % FOR each component where the stochastic potential is non null
    for m = FFTPhysics2D.StochasticPotential_function_Index{n}
        if (iscell(FFTPhysics2D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics2D.StochasticProcess_function)
                FFTPhysics2D.StochasticProcess{m_noise} = FFTPhysics2D.StochasticProcess_function{m_noise}(abs(Time1),FFTGeometry2D.X,FFTGeometry2D.Y)-FFTPhysics2D.StochasticProcess_function{m_noise}(abs(Time0),FFTGeometry2D.X,FFTGeometry2D.Y);
            end
        else
            FFTPhysics2D.StochasticProcess = FFTPhysics2D.StochasticProcess_function(abs(Time1),FFTGeometry2D.X,FFTGeometry2D.Y)-FFTPhysics2D.StochasticProcess_function(abs(Time0),FFTGeometry2D.X,FFTGeometry2D.Y);
        end
        FFTPhysics2D.StochasticPotential{n,m} = -1i*FFTPhysics2D.StochasticPotential_function{n,m}(FFTPhysics2D.StochasticProcess,FFTGeometry2D.X,FFTGeometry2D.Y); % Computing and storing the stochastic potential
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
        MatcPot{n,m} = MatcPot{n,m} + (Time1-Time0)*FFTPhysics2D.Potential{n,m}; % Computing and adding the potential operator
        MatcPot{n,m} = MatcPot{n,m} + (Time1-Time0)*FFTPhysics2D.Beta*FFTPhysics2D.Nonlinearity{n,m}; % Computing and adding the nonlinear operator
        MatcPot{n,m} = MatcPot{n,m} + (Time1-Time0)*FFTPhysics2D.Beta*FFTPhysics2D.FFTNonlinearity{n,m}; % Computing and adding the non-local nonlinear operator
        MatcPot{n,m} = MatcPot{n,m} + FFTPhysics2D.TimePotential{n,m}; % Computing and adding the time-dependent potential operator
        MatcPot{n,m} = MatcPot{n,m} + FFTPhysics2D.StochasticPotential{n,m}; % Computing and adding the stochastic potential operator
    end
end
end
function [MatcDisp] = Update_MatDispx

%% Computing the explicit nonlinearity, non-local nonlinearity and the time-depend potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the time-dependent dispersion is non null
    for m = FFTPhysics2D.TimeDispersion_function_Index{n}
        % IF the integrated time dispersion is not defined
        if (isempty(FFTPhysics2D.IntegratedTimeDispersion_function{n,m}))
            FFTPhysics2D.TimeDispersion{n,m} = ((Time1-Time0)/6)*(FFTPhysics2D.TimeDispersion_function{n,m}(abs(Time0),-1i*FFTOperators2D.Gx,0) + 4*FFTPhysics2D.TimeDispersion_function{n,m}((abs(Time0 + Time1))/2,-1i*FFTOperators2D.Gx,0) + FFTPhysics2D.TimeDispersion_function{n,m}(abs(Time1),-1i*FFTOperators2D.Gx,0)); % Computing and storing the coupled time-dependent dispersion between components
        % ELSE if the integrated time dispersion is defined
        else
            FFTPhysics2D.TimeDispersion{n,m} = -1i*(FFTPhysics2D.IntegratedTimeDispersion_function{n,m}(abs(Time1),-1i*FFTOperators2D.Gx,0) - FFTPhysics2D.IntegratedTimeDispersion_function{n,m}(abs(Time0),-1i*FFTOperators2D.Gx,0));
        end
    end
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics2D.StochasticDispersion_function_Index{n}
        if (iscell(FFTPhysics2D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics2D.StochasticProcess_function)
                FFTPhysics2D.StochasticProcess{m_noise} = FFTPhysics2D.StochasticProcess_function{m_noise}(abs(Time1),FFTGeometry2D.X,FFTGeometry2D.Y)-FFTPhysics2D.StochasticProcess_function{m_noise}(abs(Time0),FFTGeometry2D.X,FFTGeometry2D.Y);
            end
        else
            FFTPhysics2D.StochasticProcess = FFTPhysics2D.StochasticProcess_function(abs(Time1),FFTGeometry2D.X,FFTGeometry2D.Y)-FFTPhysics2D.StochasticProcess_function(abs(Time0),FFTGeometry2D.X,FFTGeometry2D.Y);
        end
        FFTPhysics2D.StochasticDispersion{n,m} = -1i*FFTPhysics2D.StochasticDispersion_function{n,m}(FFTPhysics2D.StochasticProcess,-1i*FFTOperators2D.Gx,0); % Computing and storing the stochastic potential
    end
end

% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        MatcDisp{n,m} = 0;
    end
end

%% Constructing the cell matrix using the dispersion operators
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics2D.StochasticDispersion_function_Index{n}
        MatcDisp{n,m} = MatcDisp{n,m} + FFTPhysics2D.StochasticDispersion{n,m}; % Computing and adding the dispersion operator
    end
     % FOR each component where the time-dependent dispersion is non null
    for m = FFTPhysics2D.TimeDispersion_function_Index{n}
        MatcDisp{n,m} = MatcDisp{n,m} +  FFTPhysics2D.TimeDispersion{n,m}; % Computing and adding the dispersion operator
    end
end

end

function [MatcDisp] = Update_MatDispy

%% Computing the explicit nonlinearity, non-local nonlinearity and the time-depend potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the time-dependent dispersion is non null
    for m = FFTPhysics2D.TimeDispersion_function_Index{n}
        % IF the integrated time dispersion is not defined
        if (isempty(FFTPhysics2D.IntegratedTimeDispersion_function{n,m}))
            FFTPhysics2D.TimeDispersion{n,m} = ((Time1-Time0)/6)*(FFTPhysics2D.TimeDispersion_function{n,m}(abs(Time0),0,-1i*FFTOperators2D.Gy) + 4*FFTPhysics2D.TimeDispersion_function{n,m}((abs(Time0 + Time1))/2,0,-1i*FFTOperators2D.Gy) + FFTPhysics2D.TimeDispersion_function{n,m}(abs(Time1),0,-1i*FFTOperators2D.Gy)); % Computing and storing the coupled time-dependent dispersion between components
        % ELSE if the integrated time dispersion is defined
        else
            FFTPhysics2D.TimeDispersion{n,m} = -1i*(FFTPhysics2D.IntegratedTimeDispersion_function{n,m}(abs(Time1),0,-1i*FFTOperators2D.Gy) - FFTPhysics2D.IntegratedTimeDispersion_function{n,m}(abs(Time0),0,-1i*FFTOperators2D.Gy));
        end
    end
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics2D.StochasticDispersion_function_Index{n}
        if (iscell(FFTPhysics2D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics2D.StochasticProcess_function)
                FFTPhysics2D.StochasticProcess{m_noise} = FFTPhysics2D.StochasticProcess_function{m_noise}(abs(Time1),FFTGeometry2D.X,FFTGeometry2D.Y)-FFTPhysics2D.StochasticProcess_function{m_noise}(abs(Time0),FFTGeometry2D.X,FFTGeometry2D.Y);
            end
        else
            FFTPhysics2D.StochasticProcess = FFTPhysics2D.StochasticProcess_function(abs(Time1),FFTGeometry2D.X,FFTGeometry2D.Y)-FFTPhysics2D.StochasticProcess_function(abs(Time0),FFTGeometry2D.X,FFTGeometry2D.Y);
        end
        FFTPhysics2D.StochasticDispersion{n,m} = -1i*FFTPhysics2D.StochasticDispersion_function{n,m}(FFTPhysics2D.StochasticProcess,0,-1i*FFTOperators2D.Gy); % Computing and storing the stochastic potential
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
    for m = FFTPhysics2D.StochasticDispersion_function_Index{n}
        MatcDisp{n,m} = MatcDisp{n,m} + FFTPhysics2D.StochasticDispersion{n,m}; % Computing and adding the potential operator
    end
    % FOR each component where the time-dependent dispersion is non null
    for m = FFTPhysics2D.TimeDispersion_function_Index{n}
        MatcDisp{n,m} = MatcDisp{n,m} +  FFTPhysics2D.TimeDispersion{n,m}; % Computing and adding the dispersion operator
    end
end

end
end