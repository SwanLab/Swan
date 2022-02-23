classdef LHSintergrator_StiffnessElastic < LHSintegrator

    properties (Access = private)
        geometry
    end

    methods (Access = public)
        
        function obj = LHSintergrator_StiffnessElastic(cParams)
            obj.init(cParams);
            obj.material = cParams.material;                        
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end

        function LHS = compute(obj)
            lhs   = obj.computeElementalLHS();
            LHS   = obj.assembleMatrix(lhs);
        end
        
    end
    
   methods (Access = protected)
        
%         function lhs = computeElementalLHS(obj)
%             % Belytschko page 226 (243)
%             % Is the gradient okay?
%             dShape = obj.computeGradient();
% %             if (obj.dim.ndim == 2)
% %                 nelem = obj.dim.nelem;
% %                 zer = zeros(1,3,nelem);
% %                 dShape = [dShape; zer];
% %             end
%             dvolu  = obj.mesh.computeDvolume(obj.quadrature);
%             ngaus  = obj.quadrature.ngaus;
%             nelem  = obj.mesh.nelem;
%             nnode  = obj.mesh.nnode;
%             nStre = obj.dim.nstre;
%             lhs = zeros(nnode,nnode,nelem);
%             for igaus = 1:ngaus
%                 dv(1,1,:) = dvolu(igaus,:);
%                 for iStre = 1:nStre
%                    for jStre = 1:nStre
%                       dNi = dShape(:,iStre,:,igaus);
%                       C   = obj.material.C(iStre,jStre,:);
%                       dNj = dShape(:,jStre,:,igaus);
%                       dNidNj = sum(dNi.*C.*dNj,[1,2]);
%                       lhs(iStre,jStre,:) = lhs(iStre,jStre,:) + dNidNj.*dv;
%                    end
%                 end
%             end
%         end
        
          function lhs = computeElementalLHS(obj)
            dShape = obj.geometry.dNdx;
            dvolu  = obj.mesh.computeDvolume(obj.quadrature);
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.dim.nelem;

            nstre = obj.dim.nstre;
            npe   = obj.dim.ndofPerElement;
            lhs = zeros(npe,npe,nelem);
            for igaus = 1:ngaus
                Bmat = obj.computeB(dShape, igaus);
                Cmat = obj.material.C(:,:,:,igaus);
                dV    = dvolu(igaus,:)';
                for istre = 1:nstre
                    for jstre = 1:nstre
                        Cij   = squeeze(Cmat(istre,jstre,:));                        
                        c(1,1,:) = Cij.*dV;
                        Bi = Bmat(istre,:,:); 
                        CB = bsxfun(@times,Bi,c);
                        Bj = permute(Bmat(jstre,:,:),[2 1 3]); 
                        t = bsxfun(@times,CB,Bj);                     
                        lhs = lhs + t;
                    end
                end

            end
        end

   end
    
   methods (Access = private)
  
       function B = computeB(obj,dN, igaus)
           ndim = obj.dim.ndim;
           switch ndim
               case 2
                   B = obj.computeBin2D(dN,igaus);
               case 3
                   B = obj.computeBin3D(dN,igaus);
           end
       end

       function B = computeBin2D(obj,dN, igaus)
           d = obj.dim;
           nstre          = d.nstre;
           nnode          = d.nnode;
           nelem          = d.nelem;
           nunkn          = d.nunkn;
           ndofPerElement = d.ndofPerElement;
           B = zeros(nstre,ndofPerElement,nelem);
           for i = 1:nnode
               j = nunkn*(i-1)+1;
               B(1,j,:)   = dN(1,i,:,igaus);
               B(2,j+1,:) = dN(2,i,:,igaus);
               B(3,j,:)   = dN(2,i,:,igaus);
               B(3,j+1,:) = dN(1,i,:,igaus);
           end
       end

       function B = computeBin3D(obj,dN,igaus)
           d    = obj.dim;
           B = zeros(d.nstre,d.ndofPerElement,d.nelem);
           for inode=1:d.nnode
               j = d.nunkn*(inode-1)+1;
               % associated to normal strains
               B(1,j,:)   = dN(1,inode,:,igaus);
               B(2,j+1,:) = dN(2,inode,:,igaus);
               B(3,j+2,:) = dN(3,inode,:,igaus);
               % associated to shear strain, gamma12
               B(4,j,:)   = dN(2,inode,:,igaus);
               B(4,j+1,:) = dN(1,inode,:,igaus);
               % associated to shear strain, gamma13
               B(5,j,:)   = dN(3,inode,:,igaus);
               B(5,j+2,:) = dN(1,inode,:,igaus);
               % associated to shear strain, gamma23
               B(6,j+1,:) = dN(3,inode,:,igaus);
               B(6,j+2,:) = dN(2,inode,:,igaus);
           end
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