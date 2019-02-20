classdef LevelSetFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(d)
            switch d.levelSetType
                case 'circle'
                    obj = LevelSetCircle(d);
                case 'circle_2'
                    obj = LevelSetCircle2(d);
                case 'circleInclusion'
                    obj = LevelSetWithCircleInclusion(d);
                case 'sphere'
                    obj = LevelSetSphere(d);
                case 'sphereInclusion'
                    obj = LevelSetWithSphereInclusion(d);
                case 'cylinder'
                    obj = LevelSetCylinder(d);
                case 'horizontal'
                    obj = LevelSetHorizontalInclusion(Input);
                case 'square'
                    obj = LevelSetSquareInclusion(d);
                case 'smoothSquare'
                    d.m = settings.widthSquare;                    
                    obj = LevelSetSmoothSquareInclusion(d);
                case 'rectangle'
                    obj = LevelSetRectangleInclusion(d);
                case 'smoothRectangle'
                    obj = LevelSetSmoothRectangleInclusion(d);
                case 'feasible'
                    obj = LevelSetFeasible(d);
                case 'rand'
                    obj = LevelSetRandom(d);
                case 'holes'                                       
                    obj = LevelSetWithSeveralHoles(d);
                case 'full'
                    obj = LevelSetFull(d);
                case 'horizontalFibers'
                    obj = LevelSetHorizontalFibers(d);
                otherwise
                    error('Invalid initial value of design variable.');
            end
            
        end
    end
    

end

