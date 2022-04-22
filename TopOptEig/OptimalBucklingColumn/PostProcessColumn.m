classdef PostProcessColumn < handle
    
    properties (Access = public)
        m
    end
    
    properties (Access = private)
        polygon
    end
    
    properties (Access = private)
        designVariable
        dim
        mesh
        scale
    end
    
    methods (Access = public)
        
        function obj = PostProcessColumn(cParams)
            obj.init(cParams)
        end
        
        function plotColumn(obj)
            obj.createPolygon();
            obj.createMesh();
            obj.plotMesh();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.dim            = cParams.dim;
            obj.mesh           = cParams.mesh;
            obj.scale          = cParams.scale;
        end

        function createPolygon(obj)
            A = obj.designVariable.getColumnArea;
            z = sqrt(A); 
            d     = obj.dim;
            scl   = obj.scale;
            coord = obj.mesh.coord;
            nnod     = d.nelem+1;
            dimFig   = 2;
            vertElem = 4;
            vertex = zeros(vertElem*d.nelem+1,dimFig);
            for iNod = 1:nnod-1
                vertex(2*iNod-1,:)   = [scl*z(iNod)/2 coord(iNod)];
                vertex(2*iNod,:) = [scl*z(iNod)/2 coord(iNod+1)];
            end
            vertex = obj.flip(vertex,vertElem,dimFig);
            vertex(end,:) = vertex(1,:);
            obj.polygon = polyshape(vertex);
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
            d = obj.dim;
            vertex(dimFig*d.nelem+1:vertElem*d.nelem,1) = - fliplr(vertex(1:dimFig*d.nelem,1)')';
            vertex(dimFig*d.nelem+1:vertElem*d.nelem,2) = fliplr(vertex(1:dimFig*d.nelem,2)')';
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