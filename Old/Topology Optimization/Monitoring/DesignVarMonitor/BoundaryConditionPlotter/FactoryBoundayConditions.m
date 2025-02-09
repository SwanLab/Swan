classdef FactoryBoundayConditions < handle
    
    methods (Access = public, Static)
        
        function plotter = create(shallPlot,dim,axes,mesh,bc)
            if shallPlot
                if ~isempty(dim)
                switch dim
                    case '2D'
                        plotter = BoundayConditionsPlotter_2D(axes,mesh,bc);
                    case '3D'
                        plotter = BoundayConditionsPlotter_3D(axes,mesh,bc);
                end
                    else

                        plotter = BoundayConditionsPlotter_Null();
                end
            else
                plotter = BoundayConditionsPlotter_Null();
            end
        end
        
    end
    
end