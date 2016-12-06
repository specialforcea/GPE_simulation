%% Print informations/outputs in the command window
%% INPUTS:
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var2d.m)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)

function Print_Info2d(Outputs, Method)
%% Printing informations/outputs in the command window
fprintf('----------------------------------------\n');
fprintf('Iteration %d on %d\n',Method.Iterations,Method.Max_iter); % Printing the number of iterations
% FOR each component
for n = 1:Method.Ncomponents
    fprintf('---Outputs of component %1.0f---------------\n',n);
    fprintf('Square at the origin: %8.14f\n',Outputs.phi_abs_0{n}(Outputs.Iterations)); % Printing the square of the wave function at the origin
    fprintf('x-radius mean square: %8.14f\n',Outputs.x_rms{n}(Outputs.Iterations)); % Printing the rms in the x direction
    fprintf('y-radius mean square: %8.14f\n',Outputs.y_rms{n}(Outputs.Iterations)); % Printing the rms in the y direction
    fprintf('Energy:  %8.14f\n',Outputs.Energy{n}(Outputs.Iterations)); % Printing the energy
    fprintf('Chemical potential:  %8.14f\n',Outputs.Chemical_potential{n}(Outputs.Iterations)); % Printing the chemical potential of the wave function
    fprintf('Angular momentum:  %8.14f\n',Outputs.Angular_momentum{n}(Outputs.Iterations)); % Printing the angular momentum of the wave function
    if strcmp(Method.Computation,'Ground')
        fprintf('Stopping criterion (stops at %e):  %e\n',Method.Stop_crit{2}*Method.Deltat,Method.EvolutionCriterion); % Printing the stopping criterion
    end
    % IF there are user defined local functions
    if (Outputs.User_compute_local)
       % FOR each user defined function
        for m = 1:Outputs.User_defined_number_local
            fprintf(strcat(Outputs.User_defined_names_local{m}, 32, ': %8.14f\n'), Outputs.User_defined_local{n,m}(Outputs.Iterations)) % Printing the user defined function
        end 
    end
    % IF the number of iterations is superior to 1
    if (Outputs.Iterations>1)
        Energy_decay = Outputs.Energy{n}(Outputs.Iterations)-Outputs.Energy{n}(Outputs.Iterations-1); % Computing the energy decay of the wave function
    % ELSE the number of iterations is 0
    else
        Energy_decay = 0; % Setting the energy decay to zero
    end
    fprintf('Energy evolution:  %e\n',Energy_decay); % Printing the energy of the wave function
end
fprintf('----------------------------------------\n');
% IF there are user defined local functions
if (Outputs.User_compute_global)
    % FOR each user defined function
    for m = 1:Outputs.User_defined_number_global
        fprintf(strcat(Outputs.User_defined_names_global{m}, 32, ': %8.14f\n'), Outputs.User_defined_global{m}(Outputs.Iterations)) % Printing the user defined function
    end 
end
fprintf('CPU time:  %8.2f\n',Method.Cputime); % Printing the CPU time