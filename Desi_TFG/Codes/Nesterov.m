classdef Nesterov < SGD

    properties (Access = private)
       alpha
       v
    end

    methods(Access = public)

        function obj = Nesterov(s)
            obj@SGD(s); 
            obj.alpha = s{4};
            obj.v     = 0;
        end 
    end
    
    methods (Access = protected)
        function [x,grad] = step(obj,x,e,grad,F,Xb,Yb)
            x_hat    = x + obj.alpha*obj.v;
            [~,grad] = F(x_hat,Xb,Yb);
            obj.v   = obj.alpha*obj.v - e*grad;
            x        = x + obj.v;     
        end     
    end 
end