classdef StiffnessMatrixComputer < handle
    
    properties (Access = public)       
        freeStiffnessMatrix
    end

    properties (Access = private)
        stiffnessMatrix
        elementalStiffnessMatrix
    end

    properties (Access = private)
        nElem
        length
        youngModulus
        inertiaMoment
        freeNodes
    end
    
    methods (Access = public)

        function obj = StiffnessMatrixComputer(cParams)
            obj.init(cParams)
        end

        function compute(obj)
            obj.computeElementalStiffnessMatrix();
            obj.computeStiffnessMatrix();
        end

       
        function Kfree = provideFreeStiffnessMatrix(obj)
            free = obj.freeNodes;
            K = obj.stiffnessMatrix;
            Kfree  = K(free,free);
        end

    end
    
    methods (Access = private)
        
        function obj = init(obj,cParams)
            obj.nElem         = cParams.nElem;
            obj.length        = cParams.length;
            obj.youngModulus  = cParams.youngModulus;
            obj.inertiaMoment = cParams.inertiaMoment;
            obj.freeNodes     = cParams.freeNodes;
        end

        function obj = computeElementalStiffnessMatrix(obj)
            L = obj.length;
            [c1,c2,c3,c4,c5] = obj.coeffsStiffness(L);
            Ke(1,1:4) = c1*[c2 c3 -c2 c3];
            Ke(2,1:4) = c1*[c3 c4 -c3 -c5];
            Ke(3,1:4) = c1*[-c2 -c3 c2 -c3];
            Ke(4,1:4) = c1*[c3 -c5 -c3 c4];
            obj.elementalStiffnessMatrix = Ke;
        end

        function [c1,c2,c3,c4,c5] = coeffsStiffness(obj,L)
            c1 = 1/(30*L);
            c2 = 36;
            c3 = 3*L;
            c4 = 4*L^2;
            c5 = L^2;
        end


        function computeStiffnessMatrix(obj)
            Ke = obj.elementalStiffnessMatrix;
            N = obj.nElem;
            K = sparse(2*N+2, 2*N+2);
            for iElem = 1: N
                iDof=[2*iElem-1; 2*iElem; 2*(iElem+1)-1; 2*(iElem+1)];
                K(iDof,iDof)=K(iDof,iDof)+ Ke;
            end
            obj.stiffnessMatrix = K;
        end
        
        
    end
    
end