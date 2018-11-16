classdef DensityCreatorByInitializer < DensityCreator
    
    properties (Access = private)
        
    end
    
    methods 
        
        function obj = DensityCreatorByInitializer(levFib,microProblem)
            
            shape   = microProblem.element.interpolation_u.shape; 
            conec   = microProblem.geometry.interpolation.T;
            xpoints = microProblem.geometry.interpolation.xpoints;
            
            nelem = size(conec,1);
            nnode = size(shape,1);
            quadr = microProblem.element.quadrature;
            ngaus = quadr.ngaus;
            
            yn = zeros(ngaus,nelem);
            
            for igaus = 1:ngaus
                for inode =  1:nnode
                    nodes = conec(:,inode);
                    ynode  = xpoints(nodes,2);
                    yn(igaus,:) = yn(igaus,:) + shape(inode,igaus)*ynode';
                end
            end
                            
            ls = DesignVaribleInitializer_orientedFiber.computeHorizontalFibersLevelSet(levFib,yn);
            
            dens  = zeros(nelem,ngaus);
            dens(ls > 0,:) = 0;
            dens(ls < 0,:) = 1;
            obj.density =  dens;
        end
    end
    
end

