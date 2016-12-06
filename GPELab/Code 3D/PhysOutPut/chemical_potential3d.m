%% Computation of the chemical potential
%% INPUT:
%%          Phi: Wave functions (cell array)
%%          Energy: Energy of the wave function (double)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          Geometry3D: Structure containing variables concerning the geometry of the problem in 3D (structure) (see Geometry3D_Var3d.m)
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure) (see Physics3D_Var3d.m)
%% OUTPUT:
%%          Chemical_potential: Chemical potential of the wave functions (cell array)
%% FUNCTIONS USED:
%%          L2_norm3d: To integrate the nonlinearities (line 26)

function Chemical_potential = chemical_potential3d(Phi, Energy, Method, Geometry3D, Physics3D)
%% Initialization
Chemical_potential = Energy;

%% Computation of the chemical potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        % IF the energy's nonlinearity is defined
        if (isempty(Physics3D.Nonlinearity_energy_function{n,m}) == 0)
            NL_energy = real(Physics3D.Beta*Physics3D.Nonlinearity_energy_function{n,m}(Phi,Geometry3D.X,Geometry3D.Y,Geometry3D.Z).*Phi{m}.*conj(Phi{n})); % Computing the energy of the nonlinearity with the energy's nonlinearity
            NL = real(Physics3D.Beta*Physics3D.Nonlinearity_function{n,n}(Phi,Geometry3D.X,Geometry3D.Y,Geometry3D.Z).*Phi{m}.*conj(Phi{n})); % Computing the energy of the nonlinearity with the nonlinearity of the equation
            Chemical_potential{n} = Chemical_potential{n} + Geometry3D.dx*Geometry3D.dy*Geometry3D.dz*sum(sum(sum(NL))) - Geometry3D.dx*Geometry3D.dy*Geometry3D.dz*sum(sum(sum(NL))); % Computing the chemical potential
        end
    end
end