classdef LHSintegrator_StiffnessColumn < LHSintegrator
    
    properties (Access = public)
        geometry
    end
    
    properties (Access = private)
    end
    
    properties (Access = protected)
        
    end
    
    methods (Access = public)
        
        function obj = LHSintegrator_StiffnessColumn(cParams)
            obj.init(cParams)  
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
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
    
    methods (Access = private)

        function createGeometry(obj)
            q   = obj.quadrature;
            int = obj.interpolation;
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            obj.geometry = g;
        end

%         function LHS = assembleMatrixColumn(obj,Ke)
%             d  = obj.dim;
%             N  = d.nelem;
%             K  = sparse(d.ndof,d.ndof);
%             Tdof = obj.globalConnec;   %%%%%% globalConnec: Tnod o Tdof?
%             for iElem = 1:N
%                 for iRow = 1:d.ndofPerElement
%                     iDof = Tdof(iElem,iRow);
%                     for iColum = 1:d.ndofPerElement
%                         Kij = Ke(iRow,iColum,iElem);
%                         jDof = Tdof(iElem,iColum);
%                         K(iDof,jDof) = K(iDof,jDof) + Kij;
%                     end
%                 end
%             end
%             LHS = K;
%         end

    end
    
    methods (Access = protected)

    end
    
end