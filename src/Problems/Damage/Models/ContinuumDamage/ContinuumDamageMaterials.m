classdef ContinuumDamageMaterials < handle

    methods (Static, Access = public)

        function Csec = obtainTensorSecant(C,d,r)
            degFun = 1 - d.computeFunction(r);
            Csec   = Expand(degFun,4).*C;
        end

        function Ctan = obtainTensorTangent(C,d,u,r)
            obj = ContinuumDamageMaterials;
            Csec = obj.obtainTensorSecant(C,d,r);
            dmgTangent = obj.obtainDamageTangentContribution(C,d,u,r);
            Ctan = Csec - dmgTangent;
        end

    end
    
    methods (Static, Access = private)

        function dContribution = obtainDamageTangentContribution(C,d,u,r)
            epsi   = SymGrad(u);
            sigBar = DDP(epsi,C);
            dDot   = d.computeDerivative(r);
            dContribution = Expand(dDot,4).*kronProd(sigBar,sigBar,[1 2 3 4]);
        end

    end 
end