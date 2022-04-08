classdef StiffnessMatrixComputer < handle
    
    properties (Access = public)       
        
    end

    properties (Access = private)
        connectivityMatrix
        stiffnessMatrix
        elementalStiffnessMatrix
    end

    properties (Access = private)
        nElem
       % length
        dim
        mesh
        Tnod
        youngModulus
        inertiaMoment
        freeNodes
    end
    
    methods (Access = public)

        function obj = StiffnessMatrixComputer(cParams)
            obj.init(cParams)
        end

        function compute(obj)
            obj.computeConnectivityMatrix();
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
            obj.mesh          = cParams.mesh;
            obj.dim           = cParams.dim;
            obj.Tnod          = cParams.Tnod;
            obj.youngModulus  = cParams.youngModulus;
            obj.inertiaMoment = cParams.inertiaMoment;
            obj.freeNodes     = cParams.freeNodes;
        end

        function computeConnectivityMatrix(obj)
            d = obj.dim;
            nel = obj.nElem;
            Td = zeros(nel,d.nNodE*d.nDofN);
            for iElem=1 : nel
                for a = 1 : d.nNodE
                    for j=1:d.nDofN
                        i = d.nDofN*(a-1) + j;
                        Td(iElem,i) = d.nDofN*(obj.Tnod(iElem,a)-1) + j;
                    end
                end
            end
            obj.connectivityMatrix = Td;
        end

        function obj = computeElementalStiffnessMatrix(obj)
            d = obj.dim;
            Edof  = d.nNodE*d.nDofN;
            Ke = zeros(Edof ,Edof ,obj.nElem);
            for iElem=1:obj.nElem
                l = obj.computeLength(iElem,d);
                [c1,c2,c3,c4,c5] = obj.coeffsStiffness(l);
                Ke(1,1:4,iElem) = c1*[c2 c3 -c2 c3];
                Ke(2,1:4,iElem) = c1*[c3 c4 -c3 -c5];
                Ke(3,1:4,iElem) = c1*[-c2 -c3 c2 -c3];
                Ke(4,1:4,iElem) = c1*[c3 -c5 -c3 c4];
            end
            obj.elementalStiffnessMatrix = Ke;
        end

        function l = computeLength(obj,iElem,d)
            msh = obj.mesh;
            if d.nDim == 1
                xA = msh(iElem,1);
                xB = msh(iElem+1,1);
                l = abs(xA-xB);
            elseif d.nDim == 2
                xA = msh(iElem,1);
                yA = msh(iElem,2);
                xB = msh(iElem+1,1);
                yB = msh(iElem+1,2);
                l = sqrt((xB-xA)^2+(yB-yA)^2);
            end
        end

        function [c1,c2,c3,c4,c5] = coeffsStiffness(obj,l)
            c1 = 1/(30*l);
            c2 = 36;
            c3 = 3*l;
            c4 = 4*l^2;
            c5 = l^2;
        end

        function computeStiffnessMatrix(obj)
            d  = obj.dim;
            Ke = obj.elementalStiffnessMatrix;
            N  = obj.nElem;
            K  = sparse(d.nDof,d.nDof);
            Tdof = obj.connectivityMatrix;
            Edof = d.nNodE*d.nDofN;
            for iElem = 1:N
                for iRow = 1:Edof
                    iDof = Tdof(iElem,iRow);
                    for iColum = 1:Edof
                        Kij = Ke(iRow,iColum,iElem);
                        jDof = Tdof(iElem,iColum);
                        K(iDof,jDof) = K(iDof,jDof) + Kij;
                    end
                end
            end
            obj.stiffnessMatrix = K;
        end

        
    end
    
end