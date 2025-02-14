classdef LHSintegrator_Bending < LHSintegrator

    properties (Access = public)
        elementalBendingMatrix
        bendingMatrix
    end
    
    properties (Access = private)
        youngModulus
        inertiaMoment
        designVariable
        freeNodes
    end
    
    properties (Access = private)
        geometry
    end
    
    methods (Access = public)
        
        function obj = LHSintegrator_Bending(cParams)
            obj.init(cParams);
            obj.initBending(cParams);
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end

        function LHS = compute(obj)
            obj.computeElementalBendingMatrix();
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
            obj.bendingMatrix = LHS;
        end

        function Bfree = provideFreeBendingMatrix(obj)
            free = obj.freeNodes;
            B = obj.bendingMatrix;
            Bfree  = B(free,free);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalBendingMatrix(obj)
            d = obj.dim;
            E = obj.youngModulus;
            I = obj.inertiaMoment;
            nElem = obj.mesh.nelem;
            Edof = d.ndofPerElement;
            Be = zeros(Edof ,Edof ,nElem);
            l = obj.computeLength();
            [c1,c2,c3,c4] = obj.coeffsBending(l,E,I);
            Be(1,1,:) = c1.*c2;
            Be(1,2,:) = c1.*c3;
            Be(1,3,:) = -c1.*c2;
            Be(1,4,:) = c1.*c3;
            Be(2,1,:) = c1.*c3;
            Be(2,2,:) = c1.*c4;
            Be(2,3,:) = -c1.*c3;
            Be(2,4,:) = c1.*c4/2;
            Be(3,1,:) = -c1.*c2;
            Be(3,2,:) = -c1.*c3;
            Be(3,3,:) = c1.*c2;
            Be(3,4,:) = -c1.*c3;
            Be(4,1,:) = c1.*c3;
            Be(4,2,:) = c1.*c4/2;
            Be(4,3,:) = -c1.*c3;
            Be(4,4,:) = c1.*c4;
            obj.elementalBendingMatrix = Be;
            lhs  = Be;
        end

        function lhs = computeElementalLHS(obj)
            Be = obj.elementalBendingMatrix;
            A = obj.designVariable.getColumnArea;
            d = obj.dim;
            nElem = obj.mesh.nelem;
            desVar = (A.^2)';
            Edof = d.ndofPerElement;
            lhs = zeros(Edof,Edof,nElem);
            for iElem = 1:nElem
                lhs(:,:,iElem) = desVar(iElem).*Be(:,:,iElem);
            end 
        end

        function l = computeLength(obj)
            g = obj.geometry;
            l = sum(g.dvolu,2);
        end

        function [c1,c2,c3,c4] = coeffsBending(obj,l,E,I)
            c1 = (E*I./l.^3)';
            c2 = 12*ones(1,length(l));
            c3 = (6*l)';
            c4 = (4*l.^2)';
        end

    end

    methods (Access = private)

        function obj = initBending(obj,cParams)
            obj.youngModulus   = cParams.youngModulus;
            obj.inertiaMoment  = cParams.inertiaMoment;
            obj.designVariable = cParams.designVariable; 
            obj.freeNodes      = cParams.freeNodes;
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

    end

end