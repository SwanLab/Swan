classdef DamagedMaterial < handle

    properties (Access = private)
        mesh
        baseMaterial
        damage
    end   

    methods (Access = public)

        function obj = DamagedMaterial(cParams)
            obj.init(cParams) 
        end
        
        function C = obtainNonDamagedTensor(obj)
            C = obj.baseMaterial;
        end

        function Csec = obtainTensorSecant(obj,r,rOld)
            degFun = obj.computeDegradationFun(r,rOld);
            Csec   = degFun.*obj.baseMaterial;
        end

        function Ctan = obtainTensorTangent(obj,u,r,rOld)
            Csec = obj.obtainTensorSecant(r);
            dmgTangent = obj.obtainDamageTangentContribution(u,r,rOld);
            Ctan = Csec - dmgTangent;
        end

        function d = getDamage(obj,r)
            d = obj.damage.computeFunction(r);
        end

        function q = getHardening(obj,r)
            q = obj.damage.computeHardening(r);
        end
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.baseMaterial = obj.createBaseMaterial(cParams);
            obj.damage       = obj.defineDamageLaw(cParams.hardening);
        end
        
        function mat = createBaseMaterial(obj,cParams)
            s.type    = 'ISOTROPIC';
            s.ndim    = obj.mesh.ndim;
            s.young   = ConstantFunction.create(cParams.E,obj.mesh);
            s.poisson = ConstantFunction.create(cParams.nu,obj.mesh);
            mat = Material.create(s);
        end

        function damage = defineDamageLaw(obj,s)
            s.mesh = obj.mesh;
            damage = DamageLaw(s);
        end


        function degFun = computeDegradationFun(obj,r,rOld)
            d = obj.damage.computeFunction(r,rOld);
            degFun = (1-d);
        end

        function dContribution = obtainDamageTangentContribution(obj,u,r,rOld)
            C = obj.baseMaterial;
            epsi = SymGrad(u);
            sigBar = DDP(epsi,C);

            dDot = obj.damage.computeDerivative(r,rOld);
            dContribution = dDot.*OP(sigBar,sigBar);
        end

    end 
end