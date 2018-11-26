classdef ElasticDim < handle
    
    properties (Abstract, Access = protected)
        nstre
    end
    
    methods (Access = public)
        
        function n = getNstre(obj)
            n = obj.nstre;
        end
        
    end
    
    methods (Access = protected)
        
        function strain = computeStrain(obj,d_u,idx)
            strain = zeros(obj.nstre,obj.nelem,obj.quadrature.ngaus);
            for igaus = 1:obj.quadrature.ngaus
                Bmat = obj.computeB(igaus);
                for istre=1:obj.nstre
                    for inode=1:obj.nnode
                        for idime=1:obj.dof.nunkn
                            ievab = obj.dof.nunkn*(inode-1)+idime;
                            B = squeeze(Bmat(istre,ievab,:));
                            u = d_u(idx(ievab,:));
                            strain(istre,:,igaus)=strain(istre,:,igaus)+(B.*u)';
                        end
                    end
                end
            end
        end
        
      function Bmat = computeBmat(obj)
            ngaus = obj.quadrature.ngaus;
            nnode = obj.interpolation_u.nnode;
            nunkn = obj.dof.nunkn;
            nelem = obj.nelem;
            
            Bmat = zeros(ngaus,obj.nstre,nnode*nunkn,nelem);
            for igaus = 1:ngaus
                Bmat(igaus,:,:,:) = obj.computeB(igaus);
            end
        end
        
    end
    
    methods (Abstract, Access = protected)
        computeB(obj)
    end
    
end

