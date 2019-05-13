classdef SettingsLevelSetCreatorFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(s)
            switch s.type
                case 'circle'
                    obj = LevelSetCircle(s);
                case 'circle_2'
                    obj = LevelSetCircle2(s);
                case 'circleInclusion'
                    obj = SettingsLevelSetCircleInclusion(s);
                case 'sphere'
                    obj = LevelSetSphere(s);
                case 'sphereInclusion'
                    obj = LevelSetWithSphereInclusion(s);
                case 'cylinder'
                    obj = LevelSetCylinder(s);
                case 'horizontal'
                    obj = LevelSetHorizontalInclusion(Input);
                case 'square'
                    obj = LevelSetSquareInclusion(s);
                case 'smoothSquare'
                    s.m = settings.widthSquare;                    
                    obj = LevelSetSmoothSquareInclusion(s);
                case 'rectangle'
                    obj = LevelSetRectangleInclusion(s);
                case 'smoothRectangle'
                    obj = LevelSetSmoothRectangleInclusion(s);
                case 'feasible'
                    obj = LevelSetFeasible(s);
                case 'rand'
                    obj = LevelSetRandom(s);
                case 'holes'                                       
                    obj = LevelSetWithSeveralHoles(s);
                case 'full'
                    obj = LevelSetFull(s);
                case 'horizontalFibers'
                    obj = LevelSetHorizontalFibers(s);
                case 'given'
                    obj = LevelSetGiven(s);
                otherwise
                    error('Invalid initial value of design variable.');
            end
            
        end
    end
    

end

