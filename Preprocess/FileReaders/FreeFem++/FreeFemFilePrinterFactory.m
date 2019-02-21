classdef FreeFemFilePrinterFactory < handle
        
    methods (Access = public, Static)
        
        function p = create(d)
            
            switch d.type
                case 'Identical'
                    p = FreeFemFileIdenticalPrinter();                    
                case 'InputChange'
                    p = FreeFemInputPrinter();
            end
            
        end
    end    
    
end