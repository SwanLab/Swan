classdef RMSProp < SGD

    properties (Access = private)
       rho
       r
    end

    methods(Access = public)

        function self = RMSProp(s)
            self@SGD(s); 
            self.rho = s{4};
            self.r     = 0;
        end
    end
    
    methods (Access = protected)
        function [x,grad] = step(self,x,e,grad,F,Xb,Yb)
            self.r = self.rho*self.r + (1-self.rho)*grad.*grad;
            x      = x - e./(self.r+10^-6).^0.5.*grad;
        end     
    end 
end