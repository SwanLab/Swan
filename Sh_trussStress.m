classdef Sh_trussStress < handle
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        physicalProblem
        barsLength
        designVariable
        interpolator
        qVal
        pVal
        constraintStress
        maxDisplacement
        s_q
    end
    
    methods (Access = public)
        
        function obj = Sh_trussStress(cParams)
            obj.init(cParams);
        end
        
        function computeFunction(obj)
            p = obj.physicalProblem;
            i = obj.interpolator;
            i.computeSectionInertia();
            E = p.material.E;
            I = i.sectionInertia;
            l = obj.barsLength;
            p.solve();
            s = p.stress;
            sNorm = obj.computeStressConstraintNorm(s);
            obj.constraintStress = sNorm;
            obj.value = sNorm - pi^2*E*I/(p.stiffness*l).^2; % L'Stiffness matrix no entenc massa com funcionarà aquí. :(
        end

        function computeGradient(obj)            
            p = obj.physicalProblem;
            u = p.displacement;
            RHS = obj.computeAdjointRHS(u);
            p.solveDisplacementAdjoint(RHS);
            obj.gradient = -p.Adjoint.*p.stiffnessDerivative;
        end
        
        
    end
    
    methods (Access = private)

        function obj = init(obj,cParams)
            obj.physicalProblem = cParams.phyProb;
            obj.interpolator    = cParams.interp;
            obj.maxDisplacement = cParams.maxDisplacement;
            obj.qVal            = cParams.qVal;
            obj.pVal            = cParams.pVal;
            obj.barsLength      = cParams.barsLength;
            obj.varN            = length(obj.designVariable)/2;
        end

        function sNorm = computeStressConstraintNorm(obj,s)
            nEl = size(s,3);
            q   = obj.qVal;
            p   = obj.pVal;
            s   = abs(s);
            s_e = zeros(nEl,1);
            for i = 1:nEl
                s_e(i) = (s(1,1,i)^q + s(2,2,i)^q + s(3,3,i)^q + ...
                    (2*s(1,2,i))^q + (2*s(1,3,i))^q + (2*s(2,3,i))^q)^(1/q);
            end
            obj.s_q = s_e;
            s_e     = s_e.^p;
            sNorm   = (1/nEl*sum(s_e))^(1/p);
        end

        function RHS = computeAdjointRHS(obj,u)

        end

    end
    
end