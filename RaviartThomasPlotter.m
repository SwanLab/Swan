classdef RaviartThomasPlotter < handle

    methods ( Access = public)
        
        function obj = RaviartThomasPlotter()
        end
        
        function plot(obj,s)
            figure()
            hold on

            sx = [0 1 0 0.5 0 0.5 0.25 0 0.75 0.75 0 0.25 0.25 0.5 0.25];
            sy = [0 0 1 0 0.5 0.5 0 0.25 0 0.25 0.75 0.75 0.25 0.25 0.5];

            sAF.fHandle = @(x) [x(1,:,:) ; x(2,:,:)];
            sAF.ndimf   = 2;
            sAF.mesh    = s.mesh;
            xFun = AnalyticalFunction(sAF);
            p1fun = xFun.project('P1');
            coord = p1fun.evaluate([sx; sy]);

            z = s.func.evaluate(coord);
            a = quiver(coord(1,:),coord(2,:),z(1,:),z(2,:),0);
            a.Color = [0 0 0];
            
            view(0,90)
            colorbar
            shading interp
            grid on
            s.mesh.plot();
        end
        
    end
    
    methods (Access = private)
        
    end
    
end