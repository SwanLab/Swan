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
            Csec   = degFun.*obj.baseMaterial;
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
            q = obj.damage.computeHardening(r);
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
            epsi = SymGrad(u);
            sigBar = DDP(epsi,C);
            dDot = obj.damage.computeDerivative(r);
            dContribution = dDot.*OP(sigBar,sigBar);
        end

    end 
end