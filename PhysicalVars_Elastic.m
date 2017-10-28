classdef PhysicalVars_Elastic < PhysicalVariables
    %PhysicalVars_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = ?Postprocess, SetAccess = protected)
        d_u
        strain
        stress
    end
    
    methods (Access = protected, Static)
        % !! NEEDS REVISION !!
        % Compute strains (e=B·u)
        function strain = computeStrain(d_u,dim,nnode,nelem,ngaus,idx,element)
            
            strain = zeros(dim.nstre,nelem,ngaus);
            for igaus = 1:ngaus
                Bmat = element.B.value(igaus);
                Bmat = Bmat{1,1};
                for istre=1:dim.nstre
                    for inode=1:nnode
                        for idime=1:dim.nunkn
                            ievab = dim.nunkn*(inode-1)+idime;
                            strain(istre,:,igaus)=strain(istre,:,igaus)+(squeeze(Bmat(istre,ievab,:)).*d_u(idx(ievab,:)))';
                        end
                    end
                end
            end
        end
        
        % Compute stresses
        function stres = computeStress(strain,C,ngaus,nstre)
            stres = zeros(size(strain));
            for igaus = 1:ngaus
                for istre=1:nstre
                    for jstre=1:nstre
                        stres(istre,:,igaus) = stres(istre,:,igaus) + squeeze(C(istre,jstre,:))'.*strain(jstre,:,igaus);
                    end
                end
            end
        end
    end
end

