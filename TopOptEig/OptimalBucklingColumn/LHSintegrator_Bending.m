classdef LHSintegrator_Bending < LHSintegrator
    
    properties (Access = public)
        elementalBendingMatrix
    end
    
    properties (Access = private)
        youngModulus
        inertiaMoment
        designVariable
    end
    
    properties (Access = private)
        geometry
    end
    
    methods (Access = public)
        
        function obj = LHSintegrator_Bending(cParams)
            obj.init(cParams);
            obj.initSpecific(cParams);
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assemblyBendingMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
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
            obj.elementalBendingMatrix = Be;
            lhs  = Be;
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

    end

    methods (Access = private)

        function obj = initSpecific(obj,cParams)
            obj.youngModulus   = cParams.youngModulus;
            obj.inertiaMoment  = cParams.inertiaMoment;
            obj.designVariable = cParams.designVariable; 
        end

        function createGeometry(obj)
            q   = obj.quadrature;
            int = obj.interpolation;
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            obj.geometry = g;
        end

        function LHS = assemblyBendingMatrix(obj,Be)
            d = obj.dim;
            x = obj.designVariable.value;
            nElem = obj.mesh.nelem;
            B  = sparse(d.ndof, d.ndof);
            Tdof = obj.computeDOFconnec();
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
            LHS = B;
        end

        function Td = computeDOFconnec(obj)
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
        end

    end

end