classdef DamagedMaterial < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        degradation
        baseMaterial
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = DamagedMaterial(cParams)
            obj.init(cParams)
            
        end
        
        function C = obtainTensor (obj,d)
            degFun = obj.computeDegradationFun(d);
            C = obj.createDegradedMaterial(degFun);
        end

        function C = obtainNonDamagedTensor (obj)
            C = obj.baseMaterial;
        end

        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.degradation = @(d) (1-d);
            obj.baseMaterial = obj.createBaseMaterial(cParams);
        end
        
        function mat = createBaseMaterial(obj,cParams)
            sIso.ndim = obj.mesh.ndim;
            sIso.young = ConstantFunction.create(cParams.E,obj.mesh);
            sIso.poisson = ConstantFunction.create(cParams.nu,obj.mesh);
            mat = Isotropic2dElasticMaterial(sIso);
        end
        function g = computeDegradationFun (obj,d)
            
            s.operation = @(xV) obj.degradation(d.evaluate(xV));
            s.ndimf = 1;
            g = DomainFunction(s);
        end
        
                
        function mat = createDegradedMaterial(obj,fun)
            mu    = obj.baseMaterial.createShear();
            kappa = obj.baseMaterial.createBulk();
            degM  = fun.*mu;
            degK  = fun.*kappa;
            s.shear = degM;
            s.bulk  = degK;
            s.ndim  = obj.mesh.ndim;
            mat = Isotropic2dElasticMaterial(s);
        end
        
    end
    
end