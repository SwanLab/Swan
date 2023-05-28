classdef Sh_Func_Divergence < handle
    
    properties (Access = private)
       designVariable  
    end
    
    methods (Access = public)
        
        function obj = Sh_Func_Divergence(cParams)
            obj.init(cParams)            
        end
        
        function [j,dj] = computeCostAndGradient(obj)
            j  = obj.computeCost();
            dj = obj.computeGradient();            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
        end

       function j = computeCost(obj)
           theta = obj.designVariable.thetavec;
           % Posar funció
       end

       function dj = computeGradient(obj)
           theta = obj.designVariable.thetavec;           
           % Posar funció
       end                  
        
    end
    
end