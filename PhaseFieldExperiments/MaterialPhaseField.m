classdef MaterialPhaseField < handle
    
    properties (Access = public)
        material
        Gc
        fc
    end
    
    properties (Access = private)
        mesh
        materialInterpolation
        E
        nu
    end
    
    properties (Access = private)
        mu
        kappa
    end
    

    methods (Access = public)
        
        function obj = MaterialPhaseField(cParams)
            obj.init(cParams)
        end
        
        function computeMatIso(obj,quad)
            obj.computeMatIsoParams(quad);
            obj.computeMaterial();
        end
        
        function computeMatInt(obj,cParams)
            obj.computeMatIntParams(cParams);
            obj.computeMaterial();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.materialInterpolation = cParams.materialInterpolation;
            obj.E = cParams.E;
            obj.nu = cParams.nu;
            obj.Gc = cParams.Gc;
            %obj.fc = cParams.fc;
        end
        
        function computeMatIsoParams(obj,quad)
            muVal = obj.computeMuFromYoungAndNu();
            kappaVal = obj.computeKappaFromYoungAndNu();
            
            obj.mu = muVal*ones(obj.mesh.nelem,quad.ngaus);
            obj.kappa = kappaVal*ones(obj.mesh.nelem,quad.ngaus);
        end
        
        function computeMatIntParams(obj,cParams)
            phi = cParams.phi;
            quad = cParams.quadrature;
            deriv = cParams.derivative;
            
            phiVal = phi.evaluate(quad.posgp);
            phiVal = permute(phiVal,[1 3 2]); 
            phiVal = squeezeParticular(phiVal,1);
            
            matInt  = obj.materialInterpolation.computeMatProp(phiVal,deriv);
            
            obj.kappa = matInt.kappa;
            obj.mu = matInt.mu;
        end
        
        function computeMaterial(obj)
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = obj.kappa;
            s.mu    = obj.mu;
            mat = Material.create(s);
            mat.compute(s); 
            
            obj.material = mat;
        end
        
    end
    
    methods (Access = private)
 
        function kappa = computeKappaFromYoungAndNu(obj)
            kappa = obj.E/(2*(1-obj.nu));
        end
        
        function mu = computeMuFromYoungAndNu(obj)
            mu = obj.E./(2*(1+obj.nu));
        end
        
    end
    
        
end