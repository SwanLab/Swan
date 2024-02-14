classdef StressFunctions < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        strn
        strs
    end
    
    properties (Access = private)
        mesh
        materialPhaseField
        u
        phi
        quadrature
    end
    
    methods (Access = public)
        
        function obj = StressFunctions(cParams)
            obj.init(cParams)         
        end
        
        function sV = evaluate(obj,xV)
            C = obj.computeC(xV);
            e = obj.computeStrain(xV);
            s = pagemtimes(C,e);
            s = squeezeParticular(s,2);
            sV = permute(s, [1 3 2]);
        end     
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.u   = cParams.u;
            obj.phi = cParams.phi;
            obj.materialPhaseField = cParams.materialPhaseField;            
        end

        function C = computeC(obj,xV) 
           mat =  obj.materialPhaseField;
           C   = mat.computeInterpolatedMaterial(obj.phi,xV);
        end
        
        function e = computeStrain(obj,u,xV)
            strFun = u.computeSymmetricGradient(xV);
            strFun.applyVoigtNotation();           
            e = strFun.fValues;
        end
 
    end
    
end