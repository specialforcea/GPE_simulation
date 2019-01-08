kicks = [0.5,1,2,3,4,5,10,20];

for mk=1:8
    for i=1:200
        
        
        phi = load(strcat('simulation_results/12072018kick_evolve_with_soc_no_speckle/phi_',num2str(mk),'_',num2str(i),'.mat'));
        phi = phi.phi_1;
        
        fp = sq(fourier_transform(phi,Nx,deltax));
           
      
        
        plot(X*a_0,sq(phi))
        xlabel('x/m');
        ylabel('squared wavefunction');
        legend('spin1','spin0');
        axis([-1.5e-4 1e-4 0 0.03])
        drawnow;
        saveas(gcf,strcat('simulation_results/12072018kick_evolve_with_soc_no_speckle/phi/',num2str(kicks(mk)),'/','avg_phi_',num2str(mk),'_',num2str(i),'.png'));
        
        plot(f/k_R*2*pi,fp)
        xlabel('momentum/k_R');
        ylabel('squared momentum component');
        legend('spin1','spin0');
        axis([-50 50 0 1])
        drawnow;
        saveas(gcf,strcat('simulation_results/12072018kick_evolve_with_soc_no_speckle/mom/',num2str(kicks(mk)),'/','avg_mom_',num2str(mk),'_',num2str(i),'.png'));
    end
end
    