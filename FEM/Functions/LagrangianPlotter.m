classdef LagrangianPlotter < handle

    methods ( Access = public)
        
        function obj = LagrangianPlotter()
        end
        
        function plot(obj,s)
            switch s.func.order
                case 'P0'
                    obj.plotP0Func(s)

                case {'P1','P2','P3'}
                    obj.plotLagrangianFunc(s)
            end
        end
        
    end
    
    methods (Access = private)
    
        function plotP0Func(~,s)
            p1DiscFun = s.func.project('P1D');
            p1DiscFun.plot();
        end

        
        function plotLagrangianFunc(~,s)
            for idim = 1:s.func.ndimf
                c = s.func.getCoord();
                if size(c,2) > 1
                    subplot(1,s.func.ndimf,idim);
                    hold on
                    x = c(:,1);
                    y = c(:,2);
                    z = s.func.fValues(:,idim);
                    %T = delaunay(x,y);
                    T = s.mesh.connec;
                    a = trisurf(T,x,y,z);

                    view(0,90)
                    colorbar
                    shading interp
                    grid on
                    title(['dim = ', num2str(idim)]);
                    a.EdgeColor = [0 0 0];
                    a.EdgeColor = [0 0 0];
                else
                    x = c(:,1);
                    z = s.func.fValues(:,idim);
                    plot(x,z);
                    xlim([0,1])
                end
            end
        end
        
    end
    
end