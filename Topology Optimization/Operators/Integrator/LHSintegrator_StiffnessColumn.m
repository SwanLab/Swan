classdef LHSintegrator_StiffnessColumn < LHSintegrator
    
    properties (Access = public)
        geometry
        stiffnessMatrix
    end
    
    properties (Access = private)
        freeNodes
    end
    
    properties (Access = protected)
        
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

        function Kfree = provideFreeStiffnessMatrix(obj)
            free = obj.freeNodes;
            K = obj.stiffnessMatrix;
            Kfree  = K(free,free);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            d = obj.dim;
            nElem = d.nelem; 
            Edof = d.ndofPerElement;
            Ke = zeros(Edof,Edof,nElem);
            l = obj.computeLength;
            for iElem=1:nElem
                length = l(iElem);
                [c1,c2,c3,c4,c5] = obj.coeffsStiffness(length);
                Ke(1,1:4,iElem) = c1*[c2 c3 -c2 c3];
                Ke(2,1:4,iElem) = c1*[c3 c4 -c3 -c5];
                Ke(3,1:4,iElem) = c1*[-c2 -c3 c2 -c3];
                Ke(4,1:4,iElem) = c1*[c3 -c5 -c3 c4];
            end
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
            c1 = 1/(30*l);
            c2 = 36;
            c3 = 3*l;
            c4 = 4*l^2;
            c5 = l^2;
        end        

    end
    
end