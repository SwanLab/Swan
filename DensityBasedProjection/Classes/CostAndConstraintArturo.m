classdef CostAndConstraintArturo < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = CostAndConstraintArturo(cParams)
            obj.init(cParams)                            
        end
        
        function [f0val,df0dx,df0dx2] = computeCostValueAndGradient(obj,nX,nY)
            f0val       = 0;
            df0dx       = zeros(nX*nY,1);
            df0dx2      = 0*df0dx;            
        end
        
        function [fval,dfdx,dfdx2] = computeConstraintValueAndGradient(obj,costE,costI,costD,v,dcsE1,dcsI1,dcsD1,dvs_dxs,cte)
            fval        = [costE; costI; costD; v];
            fval(1:3)   = fval(1:3)/cte;
            dfdx        = [ dcsE1(:)'; dcsI1(:)'; dcsD1(:)'; dvs_dxs(:)'];
            dfdx(1:3,:) = dfdx(1:3,:)/cte;
            dfdx2       = 0*dfdx;
            
        end
        
   
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end