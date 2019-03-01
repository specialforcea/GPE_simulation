ks = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.0,2.2];
dts = (2+ks).^2 - ks.^2 - 4;




% SOC_dressed_state;
% savepath = 'simulation_results/02212019kick_evolve_with_soc';
% kick_evolve_with_soc;


% OmegaR = 0.5;
% 
% SOC_dressed_state;
% savepath = 'simulation_results/02152019kick_evolve_with_soc/05';
% kick_evolve_with_soc;

OmegaR = 1.0;
wavepath = 'simulation_results/02132019SOC_dressed_state_8TF13/phi_1_50ms_';
savepath = 'simulation_results/03012019kick_evolve_with_soc/10';
kick_evolve_with_soc;

OmegaR = 0.5;
wavepath = 'simulation_results/02152019SOC_dressed_state_8TF13/phi_1_50ms_';
savepath = 'simulation_results/03012019kick_evolve_with_soc/05';
kick_evolve_with_soc;

OmegaR = 0.7;
wavepath = 'simulation_results/02152019SOC_dressed_state_8TF13/phi_1_50ms_';
savepath = 'simulation_results/03012019kick_evolve_with_soc/07';
kick_evolve_with_soc;







