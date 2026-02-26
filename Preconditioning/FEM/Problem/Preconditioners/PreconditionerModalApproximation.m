classdef PreconditionerModalApproximation < handle
    
    properties (Access = private)
       eigenModes    
       LHSreduced       
    end
    
    properties (Access = private)
        LHS
        nBasis
    end
    
    methods (Access = public)
        
        function obj = PreconditionerModalApproximation(cParams)
            obj.init(cParams);
            obj.computeEigenModes();
            obj.computeReducedLHS();            
        end
        
        function z = apply(obj,r)
            lhsR = obj.LHSreduced;
            phi  = obj.eigenModes;          
            rR   = phi'*r;
            zR   = lhsR\rR;
            z    = phi*zR;            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.LHS = cParams.LHS; 
            obj.nBasis = cParams.nBasis;
        end
        
        function computeEigenModes(obj)
            [phi,D]=eigs(obj.LHS,obj.nBasis,'smallestabs');
            obj.eigenModes = phi;
        end        
        
        function computeReducedLHS(obj)
            phi  = obj.eigenModes;                      
            lhs = phi'*obj.LHS*phi;           
            obj.LHSreduced = lhs;
        end
        
    end
    
end