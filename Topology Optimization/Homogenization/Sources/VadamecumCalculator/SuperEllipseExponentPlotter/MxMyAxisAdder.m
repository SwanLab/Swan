classdef MxMyAxisAdder < handle
       
    methods (Access = public, Static)
        
        function add()
            xAxis = '$m_1$';
            yAxis = '$m_2$';
            xlabel(xAxis,'Interpreter','latex');
            ylabel(yAxis,'Interpreter','latex');  
        end
        
    end

    
end