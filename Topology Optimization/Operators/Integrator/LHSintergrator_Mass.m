classdef LHSintergrator_Mass < LHSintegrator
    
    methods (Access = public)
        
        function obj = LHSintergrator_Mass(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createInterpolation();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end
        
    end
    
    methods (Access = protected)
        
       function createQuadrature(obj)
           quad = Quadrature.set(obj.mesh.type);
           quad.computeQuadrature('QUADRATIC');            
           obj.quadrature = quad;
       end        
        
        function lhs = computeElementalLHS(obj)
            shapes = obj.interpolation.shape;
            quad   = obj.quadrature;
            dvolu  = obj.mesh.computeDvolume(quad);
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.mesh.nelem;
            nnode  = obj.mesh.nnode;
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
        
    end
    
end