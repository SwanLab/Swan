classdef LHSintergrator_StiffnessElastic < LHSintegrator

    methods (Access = public)
        
        function obj = LHSintergrator_StiffnessElastic(cParams)
            obj.init(cParams)
            obj.createQuadrature();
            obj.createInterpolation();
        end

        function LHS = compute(obj)
            obj.C = obj.computeElasticityMatrix();
            lhs   = obj.computeElementalLHS();
            LHS   = obj.assembleMatrix(lhs);
        end
        
    end
    
   methods (Access = protected)
        
        function lhs = computeElementalLHS(obj)
            dShape2 = obj.interpolation.deriv;
            dShape = obj.computeGradient();
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

       function C = computeElasticityMatrix(obj)
           nstre = obj.dim.nstre;
           ngaus = obj.dim.ngaus;
           ntot  = obj.dim.nt;
           Cmat    = obj.material.C;
           CmatTot = sparse(ntot,ntot);
           dvol = obj.geometry.dvolu;
           for istre = 1:nstre
               for jstre = 1:nstre
                   for igaus = 1:ngaus
                       posI = (istre)+(nstre)*(igaus-1) : ngaus*nstre : ntot;
                       posJ = (jstre)+(nstre)*(igaus-1) : ngaus*nstre : ntot;
                       
                       Ct = squeeze(Cmat(istre,jstre,:,igaus)).*dvol(:,igaus);
                       CmatTot = CmatTot + sparse(posI,posJ,Ct,ntot,ntot);
                   end
               end
           end
       end
       
       function grad = computeGradient(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            m.type = obj.mesh.type;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(obj.quadrature.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(obj.quadrature,int);
            grad = g.cartd;
       end
       
   end
    
end