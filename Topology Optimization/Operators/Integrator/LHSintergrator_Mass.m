classdef LHSintergrator_Mass < LHSintegrator
    
    methods (Access = public)
        
        function obj = LHSintergrator_Mass(cParams)
            obj.init(cParams)
            obj.createQuadrature();
            obj.createInterpolation();
        end
        
    end
    
    methods (Access = protected)
        
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