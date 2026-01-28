classdef PlotterLatticeCross < handle

    properties (Access = private)
        patchHandle
    end

    properties (Access = private)
        mesh
        discMesh
        designVariable
        nSubdomains
        Nr
        Ntheta
        xmax
        xmin
        ymax
        ymin
        x0
        y0
        Lx
        Ly
        length
        tFrame
    end

    methods (Access = public)

        function obj = PlotterLatticeCross(cParams)
            obj.init(cParams);
            obj.createFigure();
        end

        function plot(obj)
%             rho     = obj.designVariable.fun;
            coord       = obj.updateDiscMesh();
%             funp0   = rho.project('P0');
%             rhoElem = squeeze(funp0.fValues);
%             set(obj.patchHandle,'FaceVertexAlphaData',rhoElem,'FaceAlpha','flat');
            set(obj.patchHandle, 'Vertices', coord);
%             set(obj.patchHandle,'Vertices', m.coord, 'Faces', m.connec, 'FaceColor', 'red', 'EdgeColor', 'black');
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.discMesh        = cParams.discMesh;
            obj.designVariable = cParams.fun;
            obj.nSubdomains    = cParams.nSubdomains;
            obj.Nr             = cParams.Nr  ;
            obj.Ntheta         = cParams.Ntheta;
            obj.xmax           = cParams.xmax;
            obj.xmin           = cParams.xmin;
            obj.ymax           = cParams.ymax;
            obj.ymin           = cParams.ymin;
            obj.Lx = obj.xmax - obj.xmin;
            obj.Ly = obj.ymax - obj.ymin;
            obj.length         = cParams.length;
            obj.tFrame         = cParams.tFrame;
        end

        function  coord = updateDiscMesh(obj)
            h = obj.designVariable.fValues;
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            Lx = obj.Lx;
            Ly = obj.Ly;
            coord = [];
            for jDom = 1:nY
                for iDom = 1:nX
                    ind = (jDom-1)*nX + iDom;
                    % refMesh = mesh_square_X_solid(1,h(ind),obj.Nr,obj.Ntheta);
                     refMesh = meshCrossLattice(obj.length/2,obj.tFrame,h(ind,1),h(ind,2),obj.Nr,obj.Ntheta,1);
                    coord0 = refMesh.coord;
                    s.coord(:,1) = coord0(:,1)+Lx*(iDom-1);
                    s.coord(:,2) = coord0(:,2)+Ly*(jDom-1);
                    %                     s.connec = refMesh.connec;
                    coord = [coord;s.coord];
                    %                     mIJ     = Mesh.create(s);
                    %                     plot(mIJ)
                    %                     hold on;
                    %                     mSbd{jDom,iDom} = mIJ;
                end
            end
%             s.coord = coord;
%             s.connec = obj.mesh.connec;
%             dmesh = Mesh.create(s);
        end

        function createFigure(obj)
            figure;
            set(gcf,'Pointer','arrow','NumberTitle','off');
            hold on
            axis off
            axis equal
            axes = gcf().Children;
%             obj.patchHandle = patch(axes,'Faces',obj.discMesh.connec,'Vertices',obj.discMesh.coord,...
%                 'EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
obj.patchHandle = patch(axes, 'Faces', obj.discMesh.connec, 'Vertices', obj.discMesh.coord, ...
    'EdgeColor', 'k', ...          % Changed from 'none' to 'k' (black)
    'LineStyle', '-', ...          % Changed from 'none' to '-' (solid line)
    'FaceColor', [1.0, 0.4, 0.4], ... % Optional: add a light gray face color
    'FaceLighting', 'none', ...
    'AmbientStrength', .75);
        end

    end

end