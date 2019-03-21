classdef VademecumPtensorComputer < VademecumVariablesComputer
    
    properties (Access = private)
       pNorm 
    end
    
    methods (Access = public)
       
        function obj = VademecumPtensorComputer(d)
            obj.init(d);
            obj.pNorm = d.pNorm;
        end        
        
    end
    
    methods (Access = protected)
        
        function computeAmplificators(obj,imx,imy)
            obj.createAmplificatorInput();
            d = obj.amplificatorInput;
            d.pNorm = obj.pNorm;
            ga = AmplificatorComponentsCalculator(d);
            ga.compute();
            obj.vademecumData.variables{imx,imy}.Ptensor = ga.Phomog; 
            obj.vademecumData.monomials = ga.monom;
        end          
        
    end

end