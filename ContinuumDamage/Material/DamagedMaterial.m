classdef DamagedMaterial < handle
    properties (Access = private)
        mesh
        damage
        baseMaterial
        qLaw
    end   

    methods (Access = public)
        function obj = DamagedMaterial(cParams)
            obj.init(cParams) 
        end
        
        function C = obtainNonDamagedTensor(obj)
            C = obj.baseMaterial;
        end

        function Csec = obtainTensorSecant(obj)
            degFun = obj.computeDegradationFun();
            Csec = obj.createDegradedMaterial(degFun);
        end

        function Ctan = obtainTensorTanget(obj,u)
            Csec = obj.obtainTensorSecant();
            dmgTangent = @(xV) obj.obtainDamageTangentContribution(u);
            
            op = @(xV) Csec.evaluate(xV) - dmgTangent.evaluate(xV);
            Ctan = DomainFunction.create(op,obj.mesh);
        end

        function updateMaterial (obj,qLaw)
            obj.qLaw = qLaw;
            obj.damage.updateParams(qLaw);
        end

        function d = getDamage(obj)
            d = obj.damage.computeDamage();
        end

        function q = getQ(obj)
            q = obj.damage.getQFun();
        end
    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.baseMaterial = obj.createBaseMaterial(cParams);
            obj.qLaw = cParams.qLaw;
            obj.damage = damageLaw(cParams.qLaw,obj.mesh);
        end
        
        function mat = createBaseMaterial(obj,cParams)
            s.type = 'ISOTROPIC';
            s.ndim = obj.mesh.ndim;
            s.young   = ConstantFunction.create(cParams.E,obj.mesh);
            s.poisson = ConstantFunction.create(cParams.nu,obj.mesh);
            mat = Material.create(s);
        end

        function g = computeDegradationFun(obj)
            d = @(xV)obj.damage.computeDamage();
            s.operation = @(xV) (1-d(xV));
            s.ndimf = 1;
            s.mesh  = obj.mesh;
            g = DomainFunction(s);
        end

        function dContribution = obtainDamageTangentContribution(obj,u)
            C = obj.baseMaterial;

            epsi = SymGrad(u);
            sigBar = DDP(epsi,C);

            op = @(xV) obj.damage.computeDamageDerivative(xV);
            
            dDot = DomainFunction.create(op,obj.mesh);
            dContribution = Expand(dDot).*OP(sigBar,sigBar);
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