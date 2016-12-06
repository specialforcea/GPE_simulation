%% Computation of the chemical potential
%% INPUT:
%%          Phi: Wave functions (cell array)
%%          Energy: Energy of the wave function (double)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 1D (structure) (see Geometry1D_Var1d.m)
%%          Physics1D: Structure containing variables concerning the physics of the problem in 1D (structure) (see Physics1D_Var1d.m)
%% OUTPUT:
%%          Chemical_potential: Chemical potential of the wave functions (cell array)
%% FUNCTIONS USED:
%%          L2_norm1d: To integrate the nonlinearities (line 26)

function Chemical_potential = Chemical_potential1d(Phi, Energy, Method, Geometry1D, Physics1D)
%% Initialization
Chemical_potential = Energy;

%% Computation of the chemical potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        % IF the energy's nonlinearity is defined
        if (isempty(Physics1D.Nonlinearity_energy_function{n,m}) == 0)
            NL_energy = real(Physics1D.Beta*Physics1D.Nonlinearity_energy_function{n,m}(Phi,Geometry1D.X).*Phi{m}.*conj(Phi{n})); % Computing the energy of the nonlinearity with the energy's nonlinearity
            NL = real(Physics1D.Beta*Physics1D.Nonlinearity_function{n,n}(Phi,Geometry1D.X).*Phi{m}.*conj(Phi{n})); % Computing the energy of the nonlinearity with the nonlinearity of the equation
            Chemical_potential{n} = Chemical_potential{n} + Geometry1D.dx*sum(NL) - Geometry1D.dx*sum(NL_energy); % Computing the chemical potential
        end
    end
end