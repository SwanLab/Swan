classdef LHSintergrator_Stiffness < LHSintegrator
    
    properties (Access = private)
        geometry
    end
    
    methods (Access = public)
        
        function obj = LHSintergrator_Stiffness(cParams)
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
%             dShape2 = obj.interpolation.deriv;
            dShape = obj.geometry.cartd;
            dvolu  = obj.mesh.computeDvolume(obj.quadrature);
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.mesh.nelem;
            nnode  = obj.mesh.nnode;
            lhs = zeros(nnode,nnode,nelem);
            for igaus = 1:ngaus
                dv(1,1,:) = dvolu(igaus,:);
                for iNode = 1:nnode
                   for jNode = 1:nnode
                      dNi = dShape(:,iNode,:,igaus);
                      dNj = dShape(:,jNode,:,igaus);
                      dNidNj = sum(dNi.*dNj,1);
                      lhs(iNode,jNode,:) = lhs(iNode,jNode,:) + dNidNj.*dv;
                   end
                end
            end
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
       
   end
    
end