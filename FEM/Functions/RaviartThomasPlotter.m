classdef RaviartThomasPlotter < handle

    methods ( Access = public)
        
        function obj = RaviartThomasPlotter()
        end
        
        function plot(obj,s)
            if s.mesh.ndim == 2
                obj.plot2D(s);
            elseif s.mesh.ndim == 3
                obj.plot3D(s)
            end
        end
        
    end
    
    methods (Access = private)

        function plot2D(obj,s)
            figure()
            hold on

            sx = [0, 1, 0, 0.5, 0, 0.5, 0.25, 0, 0.75, 0.75, 0, 0.25, ...
                  0.25, 0.5, 0.25, 0.125, 0, 0.875, 0.875, 0, 0.125, ...
                  0.375, 0.625, 0.375, 0.5, 0, 0, 0.125, 0.125, 0.625, ...
                  0.375, 0.5, 0.375, 0.125, 0.25, 0.125, 0.75, 0.625, ...
                  0.625, 0.125, 0.125, 0.25, 0.375, 0.25, 0.375];

            sy = [0, 0, 1, 0, 0.5, 0.5, 0, 0.25, 0, 0.25, 0.75, 0.75, ...
                  0.25, 0.25, 0.5, 0, 0.125, 0, 0.125, 0.875, 0.875, ...
                  0, 0, 0.125, 0.125, 0.375, 0.625, 0.375, 0.5, 0.375, ...
                  0.625, 0.375, 0.5, 0.125, 0.125, 0.25, 0.125, 0.125, ...
                  0.25, 0.75, 0.625, 0.625, 0.25, 0.375, 0.375];

            coord = s.mesh.computeXgauss([sx;sy]);

            z = s.func.evaluate([sx;sy]);
            s.mesh.plot();
            a = quiver(coord(1,:,:),coord(2,:,:),z(1,:,:),z(2,:,:),0.7);
            a.Color = [0 0 0];
            
            view(0,90)
            grid on
        end


        function plot3D(obj,s)
            figure()
            hold on

            sx = [0, 1, 0, 0];
            sy = [0, 0, 1, 0];
            sz = [0, 0, 0, 1];

            coord = s.mesh.computeXgauss([sx;sy;sz]);

            z = s.func.evaluate([sx;sy;sz]);
            % s.mesh.plot();
            a = quiver3(coord(1,:,:),coord(2,:,:),coord(3,:,:),z(1,:,:),z(2,:,:),z(3,:,:),1);
            a.Color = [0 0 0];
            
            view(0,90)
            grid on
        end
        
    end
    
end