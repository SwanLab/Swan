classdef FactoryBoundayConditions < handle
    
    methods (Access = public, Static)
        
        function plotter = create(shallPlot,dim,axes,mesh)
            if shallPlot
                switch dim
                    case '2D'
                        plotter = BoundayConditionsPlotter_2D(axes,mesh);
                    case '3D'
                        plotter = BoundayConditionsPlotter_3D(axes,mesh);
                end
            else
                plotter = BoundayConditionsPlotter_Null(axes,mesh);
            end
        end
        
    end
    
end