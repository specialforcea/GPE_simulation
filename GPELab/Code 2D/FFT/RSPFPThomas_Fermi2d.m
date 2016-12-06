%% Computing the Full ThomasFermi preconditioner for the Relaxation SP scheme
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          FFTPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FFT context (structure) (see FFTPhysics2D_Var2d.m)
%%          FFTGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%% OUTPUT:
%%          FPTHomasFermi: the full ThomasFermi preconditioner (structure)

function FPTHomasFermi = RSPFPThomas_Fermi2d(Psi, FFTPsi, Method, FFTPhysics2D, FFTGeometry2D)
%% Initialization of the potential and nonlinear matricial operators 
PotNLMat = cell(Method.Ncomponents);
for n = 1:Method.Ncomponents
    for m = 1:Method.Ncomponents
        PotNLMat{n,m} = (n == m)*1/Method.Deltat + FFTPhysics2D.Potential{n,m}/2 + FFTPhysics2D.TimePotentialImp{n,m}/2 + FFTPhysics2D.StochasticPotential{n,m}/2 + FFTPhysics2D.Beta*Psi{n,m}/2 + FFTPhysics2D.Beta*FFTPsi{n,m}/2;
    end
end
%% Computation of the full ThomasFermi preconditioner
% IF there are more than one component
if (Method.Ncomponents ~= 1)
    LocalFThomasFermi{1,1} = 1./PotNLMat{1,1}; % Initialization of the local full ThomasFermi preconditioner
    for n = 1:Method.Ncomponents-1
        %% Storing submatrices
        for k = 1:n
            B{k} = PotNLMat{k,n+1};
            C{k} = PotNLMat{n+1,k};
        end
        %% Computing intermediary matrices
        E = PotNLMat{n+1,n+1};
        for k = 1:n
            G{k} = zeros(FFTGeometry2D.Ny,FFTGeometry2D.Nx);
        end
        for k = 1:n
            for l = 1:n
                G{k} = G{k} + C{l}.*LocalFThomasFermi{l,k};
            end
        end
        for k = 1:n
            E = E - G{k}.*B{k};
        end
        E = 1./E;
        for k = 1:n
            F{k} = zeros(FFTGeometry2D.Ny,FFTGeometry2D.Nx);
        end
        for k = 1:n
            for l = 1:n
                F{k} = F{k} + LocalFThomasFermi{k,l}.*B{l};
            end
        end
        for k = 1:n
            for l = 1:n
                H{k,l} = E.*F{k}.*G{l};
            end
        end
        %% Computing the intermediary full Thomas Fermi preconditioner
        LocalFThomasFermi{n+1,n+1} = E;
        for k = 1:n
            LocalFThomasFermi{k,n+1} = -E.*G{k};
            LocalFThomasFermi{n+1,k} = -E.*F{k};
            for l = 1:n
                LocalFThomasFermi{k,l} = LocalFThomasFermi{k,l}+H{k,l};
            end
        end
    end
FPTHomasFermi = LocalFThomasFermi; % Storing the full Thomas Fermi preconditioner
%ELSEIF there is only one component
elseif (Method.Ncomponents == 1)
    FPTHomasFermi{1} = 1./PotNLMat{1,1}; % Storing the full Thomas Fermi preconditioner
end
        