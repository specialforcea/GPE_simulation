%% Creation of the method structure
%% INPUTS(OPTIONAL):
%%          Computation: Type of computation (string, Default:'Ground')
%%          (Must either be: 'Ground' for a groundstate computation or 'Dynamic' for a dynamic computation)
%%          Ncomponents: Number of components (double, Default: 1)
%%          Type: Type of the scheme used (string, Default: 'BESP')
%%          (Must either be: 'BEFD' for the Backward Euler Finite Difference scheme, 'BESP' for the Backward Euler SPectral scheme or 'CNSP' for the Crank-Nichoslon SPectral scheme)
%%          Deltat: Time step (double) (Default: 1e-3)
%%          Stop_time: Stopping criterion concerning the total time of computation for the dynamic computation (double, Default: 1)
%%          Stop_crit: Stopping criterion concerning the evolution of the wave functions for the groundstate computation (double, Default: 1e-8)
%%          Max_iter: Stopping criterion concerning the number of iterations for the groundstate computation (double, Default: 1e6)
%%          Preconditioner: Preconditioner used in the iterative method (string, Default: 'ThomasFermi') 
%%          (Must either be: 'Laplace' for the Laplace precondioner,'ThomasFermi'for the Thomas-Fermi preconditioner or 'None' for no preconditioner)
%%          Output: Variable containing the choice to compute outputs or no (double, Default: 1)
%%          (Must either be: 0 to not compute outputs or 1 to compute outputs)
%%          Splitting: Variable containing the choice of the type of splitting scheme (string, Default: 'Strang')
%%          (Must either be: 'Lie' for the Lie splitting scheme, 'Strang' for the Strang splitting scheme or 'Fourth' for the fourth-order splitting scheme)
%%          BESP: Variable containing the choice to use the full or the partial BESP method (double, Default: 1)
%%          (Must either be: 0 to use the full BESP method or 1 to use the partial BESP method)
%%          BEFD: Variable containing the choice to inverse the system directly or to use an iterative method in the BEFD scheme (double, Default: 1) 
%%          (Must either be: 0 to inverse the system directly or 1 to use an iterative method)
%%          Iterative_tol: Stopping criterion for the iterative method concerning the evolution of the wave functions (double, Default: 1e-9)
%%          Iterative_maxit: Stopping criterion for the iterative method concerning the number of iterations (double, Default: 1e3)
%%          Iterative_restart: Number of iterations before restarting the GMRES iterative method for the BEFD method (double, Default: [])
%% OUTPUT:
%%          Method: Structure containing variables concerning the method (structure)

function [Method] = Method_Var2d(varargin)
%% Analysis of the inputs
valid_Computation = {'Ground','Dynamic'}; % List of valid inputs for the computation type
valid_Precond = {'None','Laplace','ThomasFermi','FThomasFermi','FLaplace'}; % List of valid inputs for the precondioner
valid_Method = {'BESP','AccBESP','CNSP','BEFD','CNFD','Relaxation','Splitting'}; % List of valid inputs for the type of scheme
valid_Splitting = {'Lie','Strang','Fourth'}; % List of valid inputs for the type of splitting
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Computation','Ground',@(x)ischar(validatestring(x,valid_Computation))); % Optional input 'Computation' with default value 'Ground'
Analyse_Var.addOptional('Ncomponents',1,@(x)isposintscalar(x)); % Optional input 'Ncomponents' with default value 1
Analyse_Var.addOptional('Type','BESP',@(x)ischar(validatestring(x,valid_Method))); % Optional input 'Type' which must be part of the valid inputs with default value 'CNSP'
Analyse_Var.addOptional('Deltat',1e-3,@(x)isposrealscalar(x)); % Optional input 'Deltat' with default value 1e-2
Analyse_Var.addOptional('Stop_time',1,@(x)isposrealscalar(x)+isempty(x)); % Optional input 'Stop_time' with default value 1
Analyse_Var.addOptional('Stop_crit',{'MaxNorm',1.e-6},@(x)iscell(x)+isempty(x)); % Optional input 'Stop_crit' with default value 1e-6
Analyse_Var.addOptional('Max_iter',1e6,@(x)isposintscalar(x)+isempty(x)); % Optional input 'Max_iter' with default value 1e6
Analyse_Var.addOptional('Preconditioner','FLaplace',@(x)ischar(validatestring(x,valid_Precond))+isempty(x)); % Optional input 'Preconditioner' which must be part of the valid inputs with default value 'ThomasFermi'
Analyse_Var.addOptional('Output',1,@(x)(x==1 || x==0)); % Optional input 'Output' which must either be 1 or 0 with default value 1
Analyse_Var.addOptional('Splitting','Strang',@(x)ischar(validatestring(x,valid_Splitting))+isempty(x)); % Optional input 'Splitting' which must either be 'Lie', 'Strang' or 'Fourth'
Analyse_Var.addOptional('BESP',0,@(x)(x==1 || x==0)+isempty(x)); % Optional input 'BESP' which must either be 1 or 0 with default value 0
Analyse_Var.addOptional('Solver_FD',0,@(x)(x==1 || x==0)+isempty(x)); % Optional input 'Solver_FD' which must either be 1 or 0 with default value 0
Analyse_Var.addOptional('Iterative_tol',1e-9,@(x)isposrealscalar(x)); % Optional input 'Iterative_tol' with default value 1e-9
Analyse_Var.addOptional('Iterative_maxit',1e3,@(x)isposintscalar(x)); % Optional input 'Iterative_maxit' with default value 1e3

%% Parsing inputs and creating the Method structure
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
% Contructing the Method structure
Method.Computation = Analyse_Var.Results.Computation; %Storing the 'Computation' input
Method.Ncomponents = Analyse_Var.Results.Ncomponents; %Storing the 'Ncomponents' input
Method.Type = Analyse_Var.Results.Type; %Storing the 'Type' input
Method.Deltat = Analyse_Var.Results.Deltat; %Storing the 'Deltat' input
Method.Stop_time = Analyse_Var.Results.Stop_time; %Storing the 'Stop_time' input
Method.Stop_crit = Analyse_Var.Results.Stop_crit; %Storing the 'Stop_crit' input
Method.Max_iter = Analyse_Var.Results.Max_iter; %Storing the 'Max_iter' input
Method.Precond = Analyse_Var.Results.Preconditioner; %Storing the 'Preconditionner' input
Method.Output = Analyse_Var.Results.Output; %Storing the 'Output' input
Method.Splitting = Analyse_Var.Results.Splitting; %Storing the 'Output' input
Method.Solver_BESP = Analyse_Var.Results.BESP; %Storing the 'BESP' input
Method.Solver_FD = Analyse_Var.Results.Solver_FD; %Storing the 'Solver_FD' input
Method.Iterative_tol = Analyse_Var.Results.Iterative_tol; %Storing the 'Iterative_tol' input
Method.Iterative_maxit = Analyse_Var.Results.Iterative_maxit; %Storing the 'Iterative_maxit' input
if (isempty(Method.Stop_crit) == 0)
    Method.EvolutionCriterion = 2*Method.Stop_crit{2}; % Initializing the global evolution variable
else
    Method.EvolutionCriterion = 1; % Initializing the global evolution variable
end
Method.Iterations = 0; % Initializing the number of iterations
Method.Normalization = 'Multi'; % Normalization with respect to all the components
Method.NParticles = ones(Method.Ncomponents,1); % Number of particles in each components