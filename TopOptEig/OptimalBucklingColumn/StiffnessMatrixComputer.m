classdef StiffnessMatrixComputer < handle
    
    properties (Access = public)       
        
    end

    properties (Access = private)
        connectivityMatrix
        stiffnessMatrix
        elementalStiffnessMatrix
    end

    properties (Access = private)
        dim
        mesh
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
            obj.assemblyStiffnessMatrix();
        end

       
        function Kfree = provideFreeStiffnessMatrix(obj)
            free = obj.freeNodes;
            K = obj.stiffnessMatrix;
            Kfree  = K(free,free);
        end

    end
    
    methods (Access = private)
        
        function obj = init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.dim           = cParams.dim;
            obj.youngModulus  = cParams.youngModulus;
            obj.inertiaMoment = cParams.inertiaMoment;
            obj.freeNodes     = cParams.freeNodes;
        end

        function computeConnectivityMatrix(obj)
            Tnod = obj.mesh.connec;
            d = obj.dim;
            nElem = d.nelem;
            Td = zeros(nElem,d.nnode*d.nstre);
            for iElem=1 : nElem
                for a = 1 : d.nnode
                    for j=1:d.nstre
                        i = d.nstre*(a-1) + j;
                        Td(iElem,i) = d.nstre*(Tnod(iElem,a)-1) + j;
                    end
                end
            end
            obj.connectivityMatrix = Td;
        end

        function obj = computeElementalStiffnessMatrix(obj)
            d = obj.dim;
            Edof  = d.nnode*d.nstre;
            nElem = d.nelem;           
            Ke = zeros(Edof,Edof,nElem);
            for iElem=1:nElem
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
            coord = obj.mesh.coord;
            if d.ndim == 1
                xA = coord(iElem,1);
                xB = coord(iElem+1,1);
                l = abs(xA-xB);
            elseif d.ndim == 2
                xA = coord(iElem,1);
                yA = coord(iElem,2);
                xB = coord(iElem+1,1);
                yB = coord(iElem+1,2);
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

        function assemblyStiffnessMatrix(obj)
            d  = obj.dim;
            Ke = obj.elementalStiffnessMatrix;
            N  = obj.mesh.nelem;
            K  = sparse(d.ndof,d.ndof);
            Tdof = obj.connectivityMatrix;
            Edof = d.nnode*d.nstre;
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