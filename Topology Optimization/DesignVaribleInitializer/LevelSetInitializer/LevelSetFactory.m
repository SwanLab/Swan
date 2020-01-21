classdef LevelSetFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(d)
            switch d.type
                case 'circle'
                    obj = LevelSetCircle(d);
                case 'circleInclusion'
                    obj = LevelSetWithCircleInclusion(d);
                case 'sphere'
                    obj = LevelSetSphere(d);
                case 'sphereInclusion'
                    obj = LevelSetWithSphereInclusion(d);
                case 'cylinder'
                    obj = LevelSetCylinder(d);
                case 'horizontal'
                    obj = LevelSetHorizontalInclusion(d);
                case {'squareInclusion'}
                    obj = LevelSetSquareInclusion(d);                    
                case 'smoothSquare'
                    obj = LevelSetSmoothSquareInclusion(d);
                case 'rectangle'
                    obj = LevelSetRectangle(d);
                case 'rectangleInclusion'
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
                case 'given'
                    obj = LevelSetGiven(d);
                case 'Vigdergauz'
                    obj = LevelSetVigdergauz(d);
                otherwise
                    error('Invalid initial value of design variable.');
            end
            
        end
    end
    

end

