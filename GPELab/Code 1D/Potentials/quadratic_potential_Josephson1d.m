function [Potential] = quadratic_potential_Josephson1d(Detuning_constant, Rabi_frequency)

Potential = cell(2);
Potential{1,1} = @(x) (1/2)*x.^2 + Detuning_constant;
Potential{1,2} = @(x) Rabi_frequency;
Potential{2,1} = @(x) Rabi_frequency;
Potential{2,2} = @(x) (1/2)*x.^2;