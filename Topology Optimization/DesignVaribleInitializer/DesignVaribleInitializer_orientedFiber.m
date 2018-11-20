classdef DesignVaribleInitializer_orientedFiber < LevelSetCreator
    
    properties (Access = private)
        dir
        RotMatrix
        alpha
        
        width
        v
        levelOfFibers
        fiberPosition
    end
    
    methods
        
        function obj = DesignVaribleInitializer_orientedFiber(input)
            obj.fiberPosition = input.yn;
            obj.levelOfFibers = input.levFib;
            obj.compute(input);
            obj.computeInitialLevelSet();
        end
    end
    
    methods (Access = protected)
        
        function computeInitialLevelSet(obj)
            obj.computeLevelSet();
            obj.computeDesignVariable();
        end
        
    end
    
    
    methods (Access = private)
        function UB = computeLaminateUpperBound(obj,xc,yc)
            UB = obj.RotMatrix(2,1)*(obj.mesh.coord(:,1)-xc) + obj.RotMatrix(2,2)*(obj.mesh.coord(:,2)-yc) - (obj.width/2 -1e-6);
        end
        
        function LB = computeLaminateLowerBound(obj,xc,yc)
            LB = obj.RotMatrix(2,1)*(obj.mesh.coord(:,1)-xc) + obj.RotMatrix(2,2)*(obj.mesh.coord(:,2)-yc) + (obj.width/2 -1e-6);
        end
        
        function isVoid = isVoid(obj,s)
            
            vect = obj.v(s);
            xc = vect(1);
            yc = vect(2);
            
            UB = obj.computeLaminateUpperBound(xc,yc);
            LB = obj.computeLaminateLowerBound(xc,yc);
            isVoid = UB < 0 & LB > 0;
            
        end
        
        
        function computeLevelSet(obj)
            m = obj.levelOfFibers;
            y = obj.fiberPosition;
            period = 1/(2^m);
            phase = period/4 - mod(period/4,0.00625);
            obj.levelSet = -sin(2*pi/period*(y-phase));
        end
        
        function computeDesignVariable(obj)
            phi = obj.levelSet;
            switch obj.optimizerName
                case {'SLERP','HAMILTON-JACOBI'}
                    obj.x = phi;
                otherwise
                    initial_holes = ceil(max(phi,0))>0;
                    obj.x = obj.ini_design_value*ones(obj.lsSize);
                    obj.x(initial_holes) = obj.hole_value;
            end
        end
        
    end
end

