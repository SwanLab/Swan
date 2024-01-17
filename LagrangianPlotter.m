classdef LagrangianPlotter < handle

    methods ( Access = public)
        
        function obj = LagrangianPlotter()
        end
        
        function plot(obj,s)
            figure()
            for idim = 1:s.func.ndimf
                subplot(1,s.func.ndimf,idim);
                hold on
                
                c = s.func.getCoord();
                x = c(:,1);
                y = c(:,2);
                z = s.func.fValues;
                T = delaunay(x,y);
                a = trisurf(T,x,y,z);
                
                view(0,90)
                colorbar
                shading interp
                grid on
                title(['dim = ', num2str(idim)]);
            end
        end
        
    end
    
    methods (Access = private)
    
%         function [connec,x,y,z] = createPlotMesh(~,cParams,ielem)
%             s.coord = cParams.mesh.coord(cParams.mesh.connec(ielem,:),:);
%             s.connec = delaunay(s.coord(:,1),s.coord(:,2));
%             mesh = Mesh(s);
%             
%             r.coord = cParams.interpolation.pos_nodes(1:cParams.mesh.nnodeElem,:);
%             r.connec = delaunay(r.coord(:,1),r.coord(:,2));
%             base_mesh = Mesh(r);
%             
%             for i = 1:1
%                 base_mesh = base_mesh.remesh(2);
%                 mesh = mesh.remesh(2);
%             end
%             
%             connec = mesh.connec;
%             x = mesh.coord(:,1);
%             y = mesh.coord(:,2);
%             z = cParams.func.evaluate(base_mesh.coord');
%             z = z(:,:,ielem)';
%         end
        
    end
    
end