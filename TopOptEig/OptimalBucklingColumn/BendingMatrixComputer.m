classdef BendingMatrixComputer < handle
    
    properties (Access = public)
        elementalBendingMatrix
    end

    properties (Access = private)
        bendingMatrix
        connectivityMatrix
    end

    properties (Access = private)
        nElem
        %length
        dim
        mesh
        Tnod
        youngModulus
        inertiaMoment
        designVariable
        freeNodes
    end
    
    methods (Access = public)
        
        function obj = BendingMatrixComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.computeConnectivityMatrix();
            obj.computeElementalBendingMatrix();
            obj.computeBendingMatrix();
        end

        function Bfree = provideFreeBendingMatrix(obj)
            free = obj.freeNodes;
            B = obj.bendingMatrix;
            Bfree  = B(free,free);
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nElem          = cParams.nElem;
            % obj.length         = cParams.length;
            obj.dim            = cParams.dim;
            obj.mesh           = cParams.mesh;
            obj.Tnod           = cParams.Tnod;
            obj.youngModulus   = cParams.youngModulus;
            obj.inertiaMoment  = cParams.inertiaMoment;
            obj.designVariable = cParams.designVariable;
            obj.freeNodes      = cParams.freeNodes;
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

         
        function Be = computeElementalBendingMatrix(obj)
            % L = obj.length;
            d = obj.dim;
            E = obj.youngModulus;
            I = obj.inertiaMoment;
            Edof = d.nNodE*d.nDofN;
            Be = zeros(Edof ,Edof ,obj.nElem);
            for iElem=1: obj.nElem
                l = obj.computeLength(iElem,d);
                [c1,c2,c3,c4] = obj.coeffsBending(l,E,I);
                Be(1,1:4,iElem) = c1*[c2 c3 -c2 c3];
                Be(2,1:4,iElem) = c1*[c3 c4 -c3 c4/2];
                Be(3,1:4,iElem) = c1*[-c2 -c3 c2 -c3];
                Be(4,1:4,iElem) = c1*[c3 c4/2 -c3 c4];
            end
            obj.elementalBendingMatrix  = Be;
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

        function [c1,c2,c3,c4] = coeffsBending(obj,l,E,I)
            c1 = E*I/l^3;
            c2 = 12;
            c3 = 6*l;
            c4 = 4*l^2;
        end
        
        function computeBendingMatrix(obj)
            d = obj.dim;
            x = obj.designVariable.value;
            N = obj.nElem;
            Be = obj.elementalBendingMatrix;           
            B  = sparse(d.nDof, d.nDof);
            Tdof = obj.connectivityMatrix;
            Edof = d.nNodE*d.nDofN;
            for iElem = 1: N
                for iRow = 1:Edof
                    iDof = Tdof(iElem,iRow);
                    for iColum = 1:Edof
                        Bij = Be(iRow,iColum,iElem);
                        jDof = Tdof(iElem,iColum);
                        B(iDof,jDof) = B(iDof,jDof) + (x(iElem)^2)*Bij;
%                         iDof=[2*iElem-1; 2*iElem; 2*(iElem+1)-1; 2*(iElem+1)];
%                         B(iDof,iDof)=B(iDof,iDof)+(x(iElem)^2)*obj.elementalBendingMatrix;
    
                    end
                end
            end
            obj.bendingMatrix = B;
        end
    end

end