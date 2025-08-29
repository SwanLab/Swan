classdef DamagedMaterial < handle

    properties (Access = private)
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

        function Csec = obtainTensorSecant(obj,r)
            degFun = obj.computeDegradationFun(r);
            Csec   = Expand(degFun,4).*obj.baseMaterial;
        end

        function Ctan = obtainTensorTangent(obj,u,r)
            Csec = obj.obtainTensorSecant(r);
            dmgTangent = obj.obtainDamageTangentContribution(u,r);
            Ctan = Csec - dmgTangent;
        end

        function d = getDamage(obj,r)
            d = obj.damage.computeFunction(r);
        end

        function q = getHardening(obj,r)
            q = obj.damage.getHardening(r);
        end
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.baseMaterial = cParams.baseMaterial;
            obj.damage       = cParams.damage;
        end

        function degFun = computeDegradationFun(obj,r)
            d = obj.damage.computeFunction(r);
            degFun = (1-d);
        end

        function dContribution = obtainDamageTangentContribution(obj,u,r)
            C = obj.baseMaterial;
            epsi   = SymGrad(u);
            sigBar = DDP(epsi,C);
            dDot = obj.damage.computeDerivative(r);
            dContribution = Expand(dDot,4).*kronProd(sigBar,sigBar,[1 2 3 4]);
        end

    end 
end