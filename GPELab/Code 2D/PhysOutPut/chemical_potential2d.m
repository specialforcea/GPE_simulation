%% Computation of the chemical potential
%% INPUT:
%%          Phi: Wave functions (cell array)
%%          Energy: Energy of the wave function (double)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)
%% OUTPUT:
%%          Chemical_potential: Chemical potential of the wave functions (cell array)
%% FUNCTIONS USED:
%%          L2_norm2d: To integrate the nonlinearities (line 26)

function Chemical_potential = chemical_potential2d(Phi, Energy, Method, Geometry2D, Physics2D)
%% Initialization
Chemical_potential = Energy;

%% Computation of the chemical potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        % IF the energy's nonlinearity is defined
        if (isempty(Physics2D.Nonlinearity_energy_function{n,m}) == 0)
            NL_energy = real(Physics2D.Beta*Physics2D.Nonlinearity_energy_function{n,m}(Phi,Geometry2D.X,Geometry2D.Y).*Phi{m}.*conj(Phi{n})); % Computing the energy of the nonlinearity with the energy's nonlinearity
            NL = real(Physics2D.Beta*Physics2D.Nonlinearity_function{n,n}(Phi,Geometry2D.X,Geometry2D.Y).*Phi{m}.*conj(Phi{n})); % Computing the energy of the nonlinearity with the nonlinearity of the equation
            Chemical_potential{n} = Chemical_potential{n} + Geometry2D.dx*Geometry2D.dy*sum(sum(NL)) - Geometry2D.dx*Geometry2D.dy*sum(sum(NL_energy)); % Computing the chemical potential
        end
    end
end