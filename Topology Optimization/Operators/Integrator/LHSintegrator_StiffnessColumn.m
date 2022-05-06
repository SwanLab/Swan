classdef LHSintegrator_StiffnessColumn < LHSintegrator
    
    properties (Access = public)
        geometry
        stiffnessMatrix
    end
    
    properties (Access = private)
        freeNodes
    end
    
    methods (Access = public)
        
        function obj = LHSintegrator_StiffnessColumn(cParams)
            obj.init(cParams) 
            obj.initStiffnessColumn(cParams);
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
            obj.stiffnessMatrix = LHS;
        end

        function [Kfree,free] = provideFreeStiffnessMatrix(obj)
            free = obj.freeNodes;
            K = obj.stiffnessMatrix;
            Kfree  = K(free,free);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            d = obj.dim;
            nElem = obj.mesh.nelem; 
            Edof = d.ndofPerElement;
            Ke = zeros(Edof,Edof,nElem);
            l = obj.computeLength();
            [c1,c2,c3,c4,c5] = obj.coeffsStiffness(l);
            Ke(1,1,:) = c1.*c2;
            Ke(1,2,:) = c1.*c3;
            Ke(1,3,:) = -c1.*c2;
            Ke(1,4,:) = c1.*c3;
            Ke(2,1,:) = c1.*c3;
            Ke(2,2,:) = c1.*c4;
            Ke(2,3,:) = -c1.*c3;
            Ke(2,4,:) = -c1.*c5;
            Ke(3,1,:) = -c1.*c2;
            Ke(3,2,:) = -c1.*c3;
            Ke(3,3,:) = c1.*c2;
            Ke(3,4,:) = -c1.*c3;
            Ke(4,1,:) = c1.*c3;
            Ke(4,2,:) = -c1.*c5;
            Ke(4,3,:) = -c1.*c3;
            Ke(4,4,:) = c1.*c4;
            lhs = Ke;
        end

    end
    
    methods (Access = private)

        function initStiffnessColumn(obj,cParams)
            obj.freeNodes = cParams.freeNodes;
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

        function l = computeLength(obj)
            g = obj.geometry;
            l = sum(g.dvolu,2);
        end

        function [c1,c2,c3,c4,c5] = coeffsStiffness(obj,l)
            c1 = (1./(30*l))';
            c2 = 36*ones(1,length(l));
            c3 = (3*l)';
            c4 = (4*l.^2)';
            c5 = (l.^2)';
        end        

    end
    
end