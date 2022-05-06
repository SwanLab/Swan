classdef LHSintegrator_Mass < LHSintegrator

    properties (Access = private)
        quadType
    end

    methods (Access = public)
        
        function obj = LHSintegrator_Mass(cParams)
            obj.init(cParams);
            obj.quadType = cParams.quadType;
            obj.createQuadrature();
            obj.createInterpolation();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end
        
    end
    
    methods (Access = protected)
        
        function lhs = computeElementalLHS(obj)
            shapes = obj.interpolation.shape;
            quad   = obj.quadrature;
            dvolu  = obj.mesh.computeDvolume(quad);
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.mesh.nelem;
            nnode  = obj.mesh.nnodeElem;
            lhs = zeros(nnode,nnode,nelem);
            for igaus = 1:ngaus
                dv(1,1,:) = dvolu(igaus,:);
                Ni = shapes(:,igaus);
                Nj = shapes(:,igaus);
                NiNj = Ni*Nj';
                Aij = bsxfun(@times,NiNj,dv);
                lhs = lhs + Aij;
            end
        end

       function createQuadrature(obj)
           quad = Quadrature.set(obj.mesh.type);
           quad.computeQuadrature(obj.quadType);
           obj.quadrature = quad;
       end
        
    end
    
end