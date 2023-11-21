classdef PostProcessColumn < handle
    
    properties (Access = public)
        m
    end
    
    properties (Access = private)
        polygon
    end
    
    properties (Access = private)
        designVariable
        sectionVariables
        mesh
        scale
        frames
        makeGIF
    end
    
    methods (Access = public)
        
        function obj = PostProcessColumn(cParams)
            obj.init(cParams)
        end
        
        function plotColumn(obj)
            obj.createPolygon();
            obj.createMesh();
            obj.plotMesh();
            obj.create3Dplot();
            switch obj.makeGIF
                case 'Y'
                    obj.createGIF();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.sectionVariables = cParams.sectionVariables;
            obj.mesh           = cParams.mesh;
            obj.scale          = cParams.scale;
            obj.makeGIF        = cParams.makeGIF;
            obj.frames         = cParams.optimizer.outputFunction.monitoring.monitorDocker.designVarMonitor.frames;
        end

        function createPolygon(obj)
            z = obj.sectionVariables.computeArea();
            scl   = obj.scale;
            coord = obj.mesh.coord;
            nnod     = obj.mesh.nelem+1;
            dimFig   = 2;
            vertElem = 4;
            vertex = zeros(vertElem*obj.mesh.nelem+1,dimFig);
            for iNod = 1:nnod-1
                vertex(2*iNod-1,:)   = [scl*z(iNod)/2 coord(iNod)];
                vertex(2*iNod,:) = [scl*z(iNod)/2 coord(iNod+1)];
            end
            vertex = obj.flip(vertex,vertElem,dimFig);
            vertex(end,:) = vertex(1,:);
            obj.polygon = polyshape(vertex);
        end
        
        function create3Dplot(obj)
            nElem = obj.mesh.nelem;
            nVar = obj.designVariable.nDesignVar;
            s.designVariableValue = obj.designVariable.value(1:nVar*nElem);
            s.coord = obj.mesh.coord;
            s.type = 'cylinderBuckling'; %'cylinderBuckling'/'holedCircle'/'rectangularColumn'/'rectangularHoleColumn'
            plt = Plot3DBucklingColumn(s);
            plt.compute();
        end

        function createGIF(obj)
            fr = obj.frames;
            nImages = length(fr);
            filename = "testAnimated.gif"; 
                for idx = 1:nImages
                    im = frame2im(fr{idx});
                    [A,map] = rgb2ind(im,256);
                    if idx == 1
                        imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",1);
                    else
                        imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",1);
                    end
                end
        end
        
        function createMesh(obj)
            pgon = obj.polygon;
            tr = triangulation(pgon);
            model = createpde;
            tnodes = tr.Points';
            telements = tr.ConnectivityList';
            geometryFromMesh(model,tnodes,telements);
            generateMesh(model,'GeometricOrder','linear','Hmax',0.1)
            coord = model.Mesh.Nodes';
            connec = model.Mesh.Elements';
            s.coord = coord;
            s.connec = connec;
            obj.m = Mesh(s);
        end

        function vertex = flip(obj,vertex,vertElem,dimFig)
            nElem = obj.mesh.nelem;            
            nodes = dimFig*nElem+1:vertElem*nElem; 
            vertex(nodes,1) = - fliplr(vertex(1:dimFig*nElem,1)')';
            vertex(dimFig*nElem+1:vertElem*nElem,2) = fliplr(vertex(1:dimFig*nElem,2)')';
        end
        
        function plotMesh(obj)
            figure
            clf
            obj.m.plot();
            grid on
            grid minor
            title('Optimal clamped-clamped column','Interpreter', 'latex','FontSize',20, 'fontweight','b');
            xlabel('A(x)','Interpreter', 'latex','fontsize',14,'fontweight','b');
            ylabel('x','Interpreter', 'latex','fontsize',14,'fontweight','b');
        end
          
    end
    
end