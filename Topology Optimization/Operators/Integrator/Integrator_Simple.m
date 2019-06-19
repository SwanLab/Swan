classdef Integrator_Simple < Integrator  
    
    methods (Access = public)
        
        function obj = Integrator_Simple(cParams)
            obj.init(cParams)
        end
        
        function A = computeLHS(obj,globalConnec,npnod)
            interpolation = Interpolation.create(obj.mesh,'LINEAR');
            quadrature = obj.computeQuadrature(obj.mesh.geometryType);
            interpolation.computeShapeDeriv(quadrature.posgp);
            geometry = Geometry(obj.mesh,'LINEAR');
            geometry.computeGeometry(quadrature,interpolation);
            nelem = obj.mesh.nelem;
            Ae = zeros(interpolation.nnode,interpolation.nnode,nelem);            
            for igaus = 1:quadrature.ngaus
                for inode = 1:interpolation.nnode
                    for jnode = 1:interpolation.nnode
                        Ae(inode,jnode,:) = squeeze(Ae(inode,jnode,:)) + quadrature.weigp(igaus)*interpolation.shape(inode,igaus)...
                            *interpolation.shape(jnode,igaus)*geometry.djacob(:,igaus);
                    end
                end
            end
            
            nunkn1 = 1;
            nunkn2 = 1;
            nnode1 = size(globalConnec,2);
            nnode2 = size(globalConnec,2);
            idx1 = globalConnec';
            idx2 = globalConnec';
           
            A = sparse(npnod,npnod);
            for i = 1:nnode1*nunkn1
                for j = 1:nnode2*nunkn2
                    a = squeeze(Ae(i,j,:));
                    A = A + sparse(idx1(i,:),idx2(j,:),a,npnod,npnod);
                end
            end                
           
        end        
        
    end
    
    methods (Access = protected)
        
%         function M = integrateLHS(obj,meshUnfitted)
%             if exist('meshUnfitted','var')
%                 obj.updateMeshes(meshUnfitted);
%             end
%             M = obj.computeLHS();
%         end
%         
        function computeIntegral()          
        end
        
        
        
    end
    
    methods (Access = private)
        
      
    end
    
end