%% Creation of the print structure
%% INPUTS(OPTIONAL):
%%          Printing: Variable containing the choice to print or not to print informations on the command window during the computation (double, Default: 1)
%%          (Must either be: 1 to print or 0 to not print informations on the command window during the computation)
%%          Evo: Number of iterations before printing informations again (double, Default: 5)
%%          Draw: Variable containing the choice to draw or not to draw wave functions during the computation (double, Default: 0)
%%          (Must either be: 1 to draw or 0 to not draw wave functions during the computation)
%% OUTPUT:
%%          Print: Structure containing variables concerning the printing and drawing of informations during the computations (structure)

function [Print] = Print_Var2d(varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Printing',1,@(x)(x==1 || x==0)); % Optional input 'Printing' with default value 1
Analyse_Var.addOptional('Evo', 5,@(x)isposintscalar(x)); % Optional input 'Evo' with default value 5
Analyse_Var.addOptional('Draw',1,@(x)(x==1 || x==0)); % Optional input 'Draw' with default value 1

%% Parsing inputs and creating the Method structure
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
% Contructing the Method structure
Print.Print = Analyse_Var.Results.Printing; %Storing the 'Printing' input
Print.Evo = ceil(Analyse_Var.Results.Evo); %Storing the 'Evo' input
Print.Draw = Analyse_Var.Results.Draw; %Storing the 'Draw' input