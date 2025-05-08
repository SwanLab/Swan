classdef ContinuumDamagePlotter < handle    
    properties (Access = private)
        data
    end
    
    methods (Access = public)
        function obj = ContinuumDamagePlotter(cParams)
           obj.init(cParams)
        end

        function plotDisplacementField (obj,H)
            obj.data.displacement.field.plot;
            colorbar;
            caxis([0 1-H]);
        end

        function plotDamagesField (obj,H)
            obj.data.damage.field.plot;
            colorbar;
            caxis([0 1-H]);
        end

        function plotSelector (x,y,text)
            switch x
                case 'max r'
                    x=  obj.data.r.maxValue;
                case 'min r'
                    x = obj.data.r.minValue;
                case 'disp'
                    x = obj.data.displacement.value;
                otherwise
                    error("ContinuumDamagePlotter.m: X axis label not recognised.");            
            end
            switch y
                case 'force'
                    y = obj.data.reaction;
                case 'max damage'
                    y = obj.data.damage.maxValue;
                case 'min damage'
                    y = obj.data.damage.minValue;
                case 'max q'
                    y = obj.data.q.maxValue;
                case 'min q'
                    y = obj.data.q.minValue;
                case 'total energy'
                    y = obj.data.totalEnergy;       
                otherwise
                    error("ContinuumDamagePlotter.m: Y axis label not recognised.");
            end
            obj.plotter(x,y,text);
        end
        
    end
    methods(Access = private)
        function init(obj,cParams)
            obj.data = cParams.data;
        end
        function plotter (x,y,text)
            figure();
            plot(x,y)
            title(text);
        end
    end
end

