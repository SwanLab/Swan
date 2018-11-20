classdef LevelSetFactory < handle
    
    properties
    end
    
    methods (Access = public, Static)
        
        function obj = create(initialCase,input)
            switch initialCase
                case 'circle'
                    obj = LevelSetCircle(input);
                case 'circleInclusion'
                    obj = LevelSetWithCircleInclusion(input);
                case 'sphere'
                    obj = LevelSetSphere(input);
                case 'sphereInclusion'
                    obj = LevelSetWithSphereInclusion(input);
                case 'horizontal'
                    obj = LevelSetHorizontalInclusion(Input);
                case 'square'
                    obj = LevelSetSquareInclusion(input);
                case 'rectangle'
                    obj = LevelSetRectangleInclusion(input);
                case 'feasible'
                    obj = LevelSetFeasible(input);
                case 'rand'
                    obj = LevelSetRandom(input);
                case 'holes'
                    obj = LevelSetWithSeveralHoles(input);
                case 'full'
                    obj = LevelSetFull(input);
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

