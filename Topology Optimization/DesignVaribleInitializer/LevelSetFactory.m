classdef LevelSetFactory < handle
    
    properties
    end
    
    methods (Access = public, Static)
        
        function obj = create(initialCase,input)
            switch initialCase
                case 'circle'
                    obj = LevelSetWithCircleInclusion(input);
                case 'sphere'
                    obj = LevelSetWithSphereInclusion(input);
                case 'horizontal'
                    obj = DesignVaribleInitializer_Horizontal(input);
                case 'square'
                    obj = DesignVaribleInitializer_Square(input);
                case 'rectangle'
                    obj = DesignVaribleInitializer_Rectangle(input);
                case 'feasible'
                    obj = DesignVaribleInitializer_Feasible(input);
                case 'rand'
                    obj = DesignVaribleInitializer_Random(input);
                case 'holes'
                    obj = DesignVaribleInitializer_Holes(input);
                case 'full'
                    obj = DesignVaribleInitializer_Full(input);
                case 'orientedFiber'
                    obj = DesignVaribleInitializer_orientedFiber(input);
                case 'smoothRectangle'
                     obj = DesignVaribleInitializerRoundedRectangle(input);
                otherwise 
                    error('Invalid initial value of design variable.');
            end
            
        end
    end
    
end

