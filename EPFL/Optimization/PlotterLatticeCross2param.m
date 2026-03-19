classdef PlotterLatticeCross2param < handle

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
        figHandle
        axesHandle
    end

    methods (Access = public)

        function obj = PlotterLatticeCross2param(cParams)
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

        function gif(obj, gifName ,nIter)
            % Update geometry based on new design variables
            newCoords = obj.updateDiscMesh();
            set(obj.patchHandle, 'Vertices', newCoords);
            drawnow;
            
                obj.writeFrame(gifName,nIter);
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
            obj.x0             = cParams.x0  ;
            obj.y0             = cParams.y0  ;
            obj.Lx = obj.xmax - obj.xmin;
            obj.Ly = obj.ymax - obj.ymin;
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
                    refMesh = meshLattice2param(1,h(ind,1),h(ind,2),h(ind,2),obj.Nr,obj.Ntheta);
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
             h = obj.designVariable.fValues;
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            Lx = obj.Lx;
            Ly = obj.Ly;
            coord = [];
            connec = [];
            for jDom = 1:nY
                for iDom = 1:nX
                    ind = (jDom-1)*nX + iDom;
                    refMesh = meshLattice2param(1,h(ind,1),h(ind,2),h(ind,2),obj.Nr,obj.Ntheta);
                    coord0 = refMesh.coord;
                    s.coord(:,1) = coord0(:,1)+Lx*(iDom-1);
                    s.coord(:,2) = coord0(:,2)+Ly*(jDom-1);
                    s.connec = refMesh.connec;
                    coord = [coord;s.coord];
                    connec = [connec;s.connec+(ind-1)*refMesh.nnodes];
                                        % mIJ     = Mesh.create(s);
                    %                     plot(mIJ)
                    %                     hold on;
                                        % mSbd{jDom,iDom} = mIJ;
                end
            end
                        s.coord = coord;
                        s.connec = connec;
                        obj.discMesh = Mesh.create(s);
            % Assign handle to property and make it invisible during generation
            obj.figHandle = figure('Color', 'w', 'Visible', 'on'); 
            obj.axesHandle = axes('Parent', obj.figHandle);
            
            axis(obj.axesHandle, 'equal');
            axis(obj.axesHandle, 'off');
            hold(obj.axesHandle, 'on');

            obj.patchHandle = patch('Parent', obj.axesHandle, ...
                'Faces', obj.discMesh.connec, ...
                'Vertices', obj.discMesh.coord, ...
                'EdgeColor', 'none', 'FaceColor', 'k');
                
            % xlim(obj.axesHandle, [obj.xmin, obj.xmax]);
            % ylim(obj.axesHandle, [obj.ymin, obj.ymax]);
        end

        function writeFrame(obj, gifName, nIter)
            gifname = [gifName,'.gif'];
            f = getframe(obj.figHandle); 
            [A, map] = rgb2ind(f.cdata, 256);
            if nIter == 0 || ~exist(gifname, 'file')
                imwrite(A, map, gifname, 'gif', 'LoopCount', inf, 'DelayTime', 0.01);
            else
                imwrite(A, map,gifname, 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
            end
        end

    end

end