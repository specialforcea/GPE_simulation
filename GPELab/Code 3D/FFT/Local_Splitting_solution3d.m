%% Computation of a step in time using the full Splitting scheme
%% INPUTS:
%%          Phi: Initial wave functions in the 3D geometry for the FFT (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          FFTGeometry3D: Structure containing variables concerning the geometry of the problem in 3D in the FFT context (structure) (see FFTGeometry3D_Var3d.m)
%%          FFTPhysics3D: Structure containing variables concerning the physics of the problem in 3D in the FFT context (structure) (see FFTPhysics3D_Var3d.m)
%%          FFTOperators3D: Structure containing the derivative FFT operators (structure) (see FFTOperators3D_Var3d.m)
%% OUTPUTS:
%%          Phi: Wave functions computated with the splitting method on a single step (cell array)
%% FUNCTIONS USED:
%%          bicgstab: To compute the wave functions by an iterative method (line 48,56 and 63)
%%          operator_Full_BESP_PLaplacian3d: To compute the implicit part of the BESP scheme with the Laplace preconditioner already applied (line 48)
%%          LinearLaplace_preconditioner3d: To apply the Laplace preconditioner to the explicit part of the BESP scheme (line 49)
%%          operator_Full_BESP_PThomasFermi3d: To compute the implicit part of the BESP scheme with the Thomas-Fermi preconditioner already applied (line 56)
%%          TF_preconditioner3d: To apply the Thomas-Fermi preconditioner to the explicit part of the BESP scheme (line 57)
%%          operator_Full_BESP3d: To compute the implicit part of the BESP scheme (line 63)

function Phi = Local_Splitting_solution3d(Phi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D)

%% Computing the explicit nonlinearity and the time-dependent potential
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
        FFTPhysics3D.TimePotential{n,m} = FFTPhysics3D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing and storing the coupled time-dependent potentials between components
    end
end

%If the computation is dynamic
if (strcmp(Method.Computation,'Dynamic'))
    Method.Deltat = -1i*Method.Deltat;
end

%% Computing the Splitting scheme on a single step time
Phi = operator_Splitting3d(Phi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D);