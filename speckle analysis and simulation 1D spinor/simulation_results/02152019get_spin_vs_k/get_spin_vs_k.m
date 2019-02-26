ks = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.0,2.2];
dts = (2+ks).^2 - ks.^2 - 4;


OmegaR = 0.7;


SOC_dressed_state;
savepath = 'simulation_results/02152019kick_evolve_with_soc/07';
kick_evolve_with_soc;


OmegaR = 0.5;

SOC_dressed_state;
savepath = 'simulation_results/02152019kick_evolve_with_soc/05';
kick_evolve_with_soc;
