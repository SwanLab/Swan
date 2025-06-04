classdef Sh_Func_Divergence < handle
    
    properties (Access = private)
        designVariable
        network
    end
    
    methods (Access = public)
        
        function obj = Sh_Func_Divergence(cParams)
            obj.init(cParams)            
        end
        
        function [j,dj] = computeCostAndGradient(obj,Xb,Yb)                        
            j  = obj.computeCost(Xb,Yb);
            dj = obj.computeGradient(Yb);            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.network        = cParams.network;
        end

       function j = computeCost(obj,Xb,Yb)
           j = obj.network.forwardprop(Xb,Yb) ;%sumar terme o funció MATLAB?
       end

       function dj = computeGradient(obj,Yb)
           dj = obj.network.backprop(Yb);%sumar terme o funció MATLAB?
       end                  
        
    end
    
end