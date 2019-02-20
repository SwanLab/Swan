classdef DensityCreatorByInitializer < handle
    
    properties (Access = protected)
        density
        levelSet
    end
    
    
    methods (Access = public)
        
        function obj = DensityCreatorByInitializer(levFib,microProblem,setting)
            yP0 = obj.computeGaussVerticalPosition(microProblem);
            setting.levFib = levFib;
            setting.yn = yP0;
            setting.initial_case = 'orientedFiber';
            setting.optimizer = '';
            dvCreator = DesignVariableCreator(setting,microProblem.mesh);
            densP0 = dvCreator.getValue();
            obj.density(:,1) =  densP0;
            obj.levelSet = dvCreator.getLevelSet();
        end
        
        function d = getDensity(obj)
            d = obj.density;
        end
        
        function ls = getLevelSet(obj)
            ls = obj.levelSet;            
        end
    end
    
    methods (Static,Access = private)
        
        function yP0 = computeGaussVerticalPosition(microProblem)
            shape   = microProblem.element.interpolation_u.shape;
            conec   = microProblem.geometry.interpolation.T;
            xpoints = microProblem.geometry.interpolation.xpoints;            
            nelem = size(conec,1);
            nnode = size(shape,1);
            quadr = microProblem.element.quadrature;
            ngaus = quadr.ngaus;            
            yP0 = zeros(ngaus,nelem);            
            for igaus = 1:ngaus
                for inode =  1:nnode
                    nodes = conec(:,inode);
                    ynode  = xpoints(nodes,2);
                    yP0(igaus,:) = yP0(igaus,:) + shape(inode,igaus)*ynode';
                end
            end            
        end
        
    end
            
end

