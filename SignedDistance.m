classdef SignedDistance < handle

    properties (Access = private)
        nNodeX
        nNodeY
        dx
        mesh
        fValues
        mapTofVal
    end

    methods (Access = public)
        function obj = SignedDistance()
            obj.createMesh();
            obj.createPhi();
            obj.createMapToValOperator();
        end

        function [phi,nx,ny,dX] = getPhi(obj)
            phi = obj.fValues;
            nx  = obj.nNodeX;
            ny  = obj.nNodeY;
            dX  = obj.dx;
        end

        function plotSignedDistance(obj,mapU)
            uVal = obj.mapTofVal(mapU);
            f = LagrangianFunction.create(obj.mesh,1,'P1');
            f.setFValues(uVal);
            plot(f);
            view(45,20);
        end
    end

    methods (Access = private)
        function createMesh(obj)
            % DATA
            nElemx = 100;
            nElemy = 100;
            Lx     = 2;
            Ly     = 2;

            obj.nNodeX = nElemx + 1;
            obj.nNodeY = nElemy + 1;
            obj.mesh = QuadMesh(Lx,Ly,nElemx,nElemy);
            obj.dx = Lx/nElemx;
        end

        function createPhi(obj)
%             % CIRCLE
%             s.type = 'Circle';
%             s.radius = 0.5;
%             s.xCoorCenter = 1;
%             s.yCoorCenter = 1;
%             g = GeometricalFunction(s);
%             ls = g.computeLevelSetFunction(obj.mesh);
%             obj.fValues = ls.fValues;
%             ls.print('LevelSet');

%             % SQUARE
%             s.type = 'Square';
%             s.length = 1;
%             s.xCoorCenter = 1;
%             s.yCoorCenter = 1;
%             g = GeometricalFunction(s);
%             ls = g.computeLevelSetFunction(obj.mesh);
%             obj.fValues = ls.fValues;
%             ls.print('LevelSet');

%             % RECTANGLE ROTATED
%             s.type = 'RectangleRotated';
%             s.xSide = 1.5;
%             s.ySide = 0.5;
%             s.xCoorCenter = 1;
%             s.yCoorCenter = 1;
%             s.omegaDeg = 45;
%             g = GeometricalFunction(s);
%             ls = g.computeLevelSetFunction(obj.mesh);
%             obj.fValues = ls.fValues;
%             ls.print('LevelSet');

             % FIBERS
            s.type = 'HorizontalNFibers';
            s.nFibers = 5;
            s.minyCoor = 0.1;
            s.maxyCoor = 0.9;
            g = GeometricalFunction(s);
            ls = g.computeLevelSetFunction(obj.mesh);
            obj.fValues = ls.fValues;
            ls.print('LevelSet');
        end

        function createMapToValOperator(obj)
            obj.mapTofVal = @(map) reshape(map,[obj.nNodeY*obj.nNodeX,1]);
        end
    end
end