classdef DesignVaribleInitializer_orientedFiber < DesignVaribleInitializer
        
        properties (Access = private)
            dir
            RotMatrix
            alpha

            width
            v
            LevelOfFibers
        end
        
        methods
        
        function obj = DesignVaribleInitializer_orientedFiber(settings,mesh,...
                        epsilon,dir,levFib)
            obj@DesignVaribleInitializer(settings,mesh,epsilon);
            obj.dir = dir;
            obj.LevelOfFibers = levFib;
        end
        
        function x = compute_initial_x(obj)            
%             obj.alpha = atan(obj.dir(2)/obj.dir(1));
%             obj.RotMatrix = [cos(obj.alpha) sin(obj.alpha); -sin(obj.alpha) cos(obj.alpha)];
%             
%             
%             LevelofFibers = obj.LevelOfFibers; 
%             s = linspace(-1,1,1+LevelofFibers*4);
%             
%             volumen = 0.5;
%             nFibers = length(s);
%             obj.width = (1-volumen)/nFibers;
%             
%             center  = [0.5;0.5];            
%             n = [-obj.dir(2),obj.dir(1)]';
%             smax = min(abs(0.5./n));
%             
%             obj.v = @(s) center + s*smax*n;
% 
%             
%             
%             isReallyVoid = false(size(obj.mesh.coord(:,2)));
%             for iFibers = 1:nFibers
%                 isVoid = obj.isVoid(s(iFibers));
%                 isReallyVoid = isReallyVoid | isVoid;
%             end
%             
%             obj.x(isReallyVoid) = obj.hole_value;
%             x = obj.x;
            
            yn = obj.mesh.coord(:,2);
            lev = obj.LevelOfFibers;
            phi = obj.computeHorizontalFibersLevelSet(lev,yn);
            obj.x = phi;
            %obj.x(phi>0) = obj.hole_value ;
            x = obj.x;
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
            
        end
        
        
        methods (Access = public,Static)
            
            function phi = computeHorizontalFibersLevelSet(levelOfFibers,y)
                m = levelOfFibers;
                period = 1/(2^m);
                phase = period/4 - mod(period/4,0.00625);
                phi = -sin(2*pi/period*(y-phase));
            end
            
        end
end

