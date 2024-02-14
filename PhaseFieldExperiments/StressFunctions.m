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
        
        function fxV = evaluate(obj)
            fxV = obj.strs.fValues;
        end     
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.materialPhaseField = cParams.materialPhaseField;
            obj.u = cParams.u;
            obj.phi = cParams.phi;
            obj.quadrature = cParams.quadrature;
            
            obj.computeStrain(obj.u,obj.quadrature);
            obj.computeStress(obj.phi,obj.quadrature);
        end
        
        function computeStrain(obj,u,quad)
            strFun = u.computeSymmetricGradient(quad);
            strFun.applyVoigtNotation();
            obj.strn = strFun;
        end
        
        function computeStress(obj,phi,quad)
            strainVal  = permute(obj.strn.fValues,[1 3 2]);
            strn2(:,1,:,:) = strainVal;
            
            obj.materialPhaseField.computeInterpolatedMaterial(phi,quad)
            stressVal =squeezeParticular(pagemtimes(obj.materialPhaseField.material.C,strn2),2);
            stressVal = permute(stressVal, [1 3 2]);
            
            z.mesh       = obj.mesh;
            z.fValues    = stressVal;
            z.quadrature = quad;
            obj.strs = FGaussDiscontinuousFunction(z);
        end
        
    end
    
end