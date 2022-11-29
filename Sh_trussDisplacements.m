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
        u_q
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
            u = p.displacement;
            RHS = obj.computeAdjointRHS(u);
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
            obj.u_q = u_i;
            u_i     = u_i.^p;
            uNorm   = (1/nNodes*sum(u_i))^(1/p);
        end

        function RHS = computeAdjointRHS(obj,u)
            gu     = obj.constraintU;
            nNodes = size(u,3);
            p      = obj.pVal;
            q      = obj.qVal;
            uq     = obj.u_q;
            k      = gu^(1-p)*1/nNodes*uq.^(p-q);
            RHS    = zeros(nNodes*3,1);
            for i = 1:length(k)
                uNode = k(i).*u(:,:,i);
                n1    = 3*i-2;
                n2    = 3*i;
                RHS(n1:n2) = uNode;
            end
        end

    end
    
end