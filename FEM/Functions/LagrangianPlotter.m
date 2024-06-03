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
                % plot(1,s.func.ndimf,idim);
                %hold on
                figure()
                c = s.func.getCoord();
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
                a.EdgeAlpha = 0;
    
            end
        end
        
    end
    
end