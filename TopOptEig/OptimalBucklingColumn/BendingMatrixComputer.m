classdef BendingMatrixComputer < handle
    
    properties (Access = public)
        elementalBendingMatrix
    end

    properties (Access = private)
        bendingMatrix
        connectivityMatrix
    end

    properties (Access = private)
        dim
        mesh
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
            obj.assemblyBendingMatrix();
        end

        function Bfree = provideFreeBendingMatrix(obj)
            free = obj.freeNodes;
            B = obj.bendingMatrix;
            Bfree  = B(free,free);
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.dim            = cParams.dim;
            obj.mesh           = cParams.mesh;
            obj.youngModulus   = cParams.youngModulus;
            obj.inertiaMoment  = cParams.inertiaMoment;
            obj.designVariable = cParams.designVariable;
            obj.freeNodes      = cParams.freeNodes;
        end

        function computeConnectivityMatrix(obj)
            connec = obj.mesh.connec;            
            d = obj.dim;
            nElem = obj.mesh.nelem;
            Td = zeros(nElem,d.nnode*d.nstre);
            for iElem= 1 : nElem
                for iNode = 1 : d.nnode
                    for j=1:d.nstre
                        i = d.nstre*(iNode-1) + j;
                        node = connec(iElem,iNode);
                        Td(iElem,i) = d.nstre*(node-1) + j;
                    end
                end
            end
            obj.connectivityMatrix = Td;
        end

         
        function Be = computeElementalBendingMatrix(obj)
            d = obj.dim;
            E = obj.youngModulus;
            I = obj.inertiaMoment;
            nElem = d.nelem;            
            Edof = d.nnode*d.nstre;
            Be = zeros(Edof ,Edof ,nElem);
            for iElem = 1:nElem
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

        function [c1,c2,c3,c4] = coeffsBending(obj,l,E,I)
            c1 = E*I/l^3;
            c2 = 12;
            c3 = 6*l;
            c4 = 4*l^2;
        end
        
        function assemblyBendingMatrix(obj)
            d = obj.dim;
            x = obj.designVariable.value;
            nElem = obj.mesh.nelem;
            Be = obj.elementalBendingMatrix;           
            B  = sparse(d.ndof, d.ndof);
            Tdof = obj.connectivityMatrix;
            Edof = d.nnode*d.nstre;
            for iElem = 1: nElem
                for iRow = 1:Edof
                    iDof = Tdof(iElem,iRow);
                    for iColum = 1:Edof
                        Bij = Be(iRow,iColum,iElem);
                        jDof = Tdof(iElem,iColum);
                        B(iDof,jDof) = B(iDof,jDof) + (x(iElem)^2)*Bij;
    
                    end
                end
            end
            obj.bendingMatrix = B;
        end
    end

end