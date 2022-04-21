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
        geometry
        youngModulus
        inertiaMoment
        designVariable
        freeNodes
    end
    
    methods (Access = public)
        
        function obj = BendingMatrixComputer(cParams)
            obj.init(cParams);
            obj.createGeometry();
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

        function createGeometry(obj)
            s.mesh = obj.mesh;

            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');

            q   = quad;
            int = Interpolation.create(obj.mesh,'LINEAR');
            int.computeShapeDeriv(q.posgp);

            obj.geometry = Geometry.create(s);
            obj.geometry.computeGeometry(q,int)
        end

        function computeConnectivityMatrix(obj)
            Tnod = obj.mesh.connec;
            d = obj.dim;
            nElem = d.nelem;
            Td = zeros(nElem,d.ndofPerElement);  
            ndofn = d.ndimField;
            for iElem=1 : nElem
                for a = 1 : ndofn
                    for j=1:ndofn
                        i = ndofn*(a-1) + j;
                        Td(iElem,i) = ndofn*(Tnod(iElem,a)-1) + j;
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
            Edof = d.ndofPerElement;
            Be = zeros(Edof ,Edof ,nElem);
            l = obj.computeLength();
            for iElem = 1:nElem
                length = l(iElem);
                [c1,c2,c3,c4] = obj.coeffsBending(length,E,I);
                Be(1,1:4,iElem) = c1*[c2 c3 -c2 c3];
                Be(2,1:4,iElem) = c1*[c3 c4 -c3 c4/2];
                Be(3,1:4,iElem) = c1*[-c2 -c3 c2 -c3];
                Be(4,1:4,iElem) = c1*[c3 c4/2 -c3 c4];
            end
            obj.elementalBendingMatrix  = Be;
        end

        function l = computeLength(obj)
            g = obj.geometry;
            l = sum(g.dvolu,2);
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
            Edof = d.ndofPerElement;
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