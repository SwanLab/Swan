classdef DamagedMaterial < handle
    properties (Access = private)
        mesh
        damage
        baseMaterial
    end   

    methods (Access = public)
        function obj = DamagedMaterial(cParams)
            obj.init(cParams) 
        end
        
        function C = obtainNonDamagedTensor(obj)
            C = obj.baseMaterial;
        end

        function Csec = obtainTensorSecant(obj,d)
            C = obj.baseMaterial;
            d = obj.damage;
            Csec = (1-d)*C;
            degFun = obj.computeDegradationFun(d);
            C = obj.createDegradedMaterial(degFun);
        end

        function C = obtainTensorTanget(obj,d)
            Csec = obj.obtainTensorSecant();
            dmgTangent = obj.obtainDamageTangentContribution();
            Ctan = Csec - dmgTangent;
        end

    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.degradation = @(d) (1-d);
            obj.baseMaterial = obj.createBaseMaterial(cParams);
        end
        
        function mat = createBaseMaterial(obj,cParams)
            s.type = 'ISOTROPIC';
            s.ndim = obj.mesh.ndim;
            s.young   = ConstantFunction.create(cParams.E,obj.mesh);
            s.poisson = ConstantFunction.create(cParams.nu,obj.mesh);
            mat = Material.create(s);
        end

        function g = computeDegradationFun(obj,d)
            s.operation = @(xV) obj.degradation(d.evaluate(xV));
            s.ndimf = 1;
            s.mesh  = obj.mesh;
            g = DomainFunction(s);
        end

        function dg = computeDerivativeDegradationFun(obj,d)

        end
                
        function mat = createDegradedMaterial(obj,fun)
            mu    = obj.baseMaterial.createShear(obj.mesh);
            kappa = obj.baseMaterial.createBulk(obj.mesh);
            degM  = fun.*mu;
            degK  = fun.*kappa;
            s.shear = degM;
            s.bulk  = degK;
            s.ndim  = obj.mesh.ndim;
            mat = Isotropic2dElasticMaterial(s);
        end     
    end 
end