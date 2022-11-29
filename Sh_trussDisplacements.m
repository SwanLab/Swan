classdef Sh_trussDisplacements < handle
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        physicalProblem
        designVariable
        normVal
        qVal
        pVal
        constraintU
        maxDisplacement
    end
    
    methods (Access = public)
        
        function obj = Sh_trussDisplacements(cParams)
            obj.init(cParams);
        end
        
        
        function computeFunction(obj)
            p = obj.physicalProblem;
            p.solve();
            u = p.displacement;
            uNorm = obj.computeDisplacementConstraintNorm(u);
            obj.constraintU = uNorm;
            obj.value = uNorm - obj.maxDisplacement;
        end

        function computeGradient(obj)            
            p = obj.physicalProblem;
            RHS = obj.computeAdjointRHS();
            p.solveDisplacementAdjoint(RHS);
            obj.gradient = -p.Adjoint.*p.stiffnessDerivative;
        end
        
        
    end
    
    methods (Access = private)

        function obj = init(obj,cParams)
            obj.physicalProblem = cParams.phyProb;
            obj.maxDisplacement = cParams.maxDisplacement;
            obj.normVal         = cParams.normVal;
            obj.qVal            = cParams.qVal;
            obj.pVal            = cParams.pVal;
            obj.varN            = length(obj.designVariable)/2;
        end

        function uNorm = computeDisplacementConstraintNorm(obj,u)
            nNodes = size(u,3);
            q      = obj.qVal;
            p      = obj.pVal;
            u_i    = zeros(1,nNodes);
            for i = 1:nNodes
                u_i(i) = (u(1,1,i)^q + u(1,2,i)^q + u(1,3,i)^q)^(1/q);
            end
            u_i = u_i.^p;
            uNorm = (1/nNodes*sum(u_i))^(1/p);
        end

        function RHS = computeAdjointRHS(obj)
            gu     = obj.constraintU;
            nNodes = size(u,3);
            p      = obj.pVal;
            RHS    = gu^(1-p)*1/nNodes; % PER ACABAR
        end

    end
    
end