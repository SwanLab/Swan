classdef Sh_Func_L2norm < handle
    
    properties (Access = private)
       designVariable  
    end
    
    methods (Access = public)
        
        function obj = Sh_Func_L2norm(cParams)
            obj.init(cParams)            
        end
        
        function [j,dj] = computeFunctionAndGradient(obj, x)
            obj.designVariable.thetavec = x;
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
           j = 0.5*(theta)*theta';
       end

       function dj = computeGradient(obj)
           theta = obj.designVariable.thetavec;           
           dj = theta;
       end                  
        
    end
    
end