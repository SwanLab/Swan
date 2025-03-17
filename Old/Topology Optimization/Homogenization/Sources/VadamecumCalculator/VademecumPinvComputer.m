classdef VademecumPinvComputer < VademecumVariablesComputer
    
    methods (Access = public)
       
        function obj = VademecumPinvComputer(d)
            obj.init(d);
        end        
        
    end
    
    methods (Access = protected)
        
        function computeAmplificators(obj,imx,imy)
            obj.createAmplificatorInput();
            d = obj.amplificatorInput;
            ga = AmplificatorSquareCalculator(d);
            ga.compute();
            obj.vademecumData.variables{imx,imy}.Ptensor = ga.Phomog;             
            obj.vademecumData.variables{imx,imy}.Pinvtensor = ga.PhomogInv();            
        end          
        
    end

end