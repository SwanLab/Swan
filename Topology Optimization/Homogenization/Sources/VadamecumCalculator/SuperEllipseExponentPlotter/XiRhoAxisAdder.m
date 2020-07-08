classdef XiRhoAxisAdder < handle
       
    methods (Access = public, Static)
        
        function add()
            xAxis = '$\xi$';
            yAxis = '$\rho$';
            xlabel(xAxis,'Interpreter','latex');
            set(gca,'xtick',[0:pi/8:pi/2]);
            set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})                         
            ylim([0 1])
            ylabel(yAxis,'Interpreter','latex');  
        end
        
    end

    
end