classdef RMSProp < SGD

    properties (Access = private)
       rho
       r
    end

    methods(Access = public)

        function obj = RMSProp(s)
            obj@SGD(s); 
            obj.rho = s{4};
            obj.r     = 0;
        end
    end
    
    methods (Access = protected)
        function [x,grad] = step(obj,x,e,grad,F,Xb,Yb)
            obj.r = obj.rho*obj.r + (1-obj.rho)*grad.*grad;
            x      = x - e./(obj.r+10^-6).^0.5.*grad;
        end     
    end 
end