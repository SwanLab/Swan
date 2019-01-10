classdef ResultsPrinterFactory < handle
        
    methods (Access = public, Static)
        function p = create(resultCase,d)
            
            switch resultCase 
                case 'Elasticity'
                    p = ElasticityResultsPrinter(d);     
                case 'TopOptProblem'
                    p = TopOptResultsPrinter(d);
            end            
            
        end
        
    end
end