classdef PerformanceTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        
    end


    methods (Test, TestTags = {'Performance'})

        function test2D(testCase)
            index = 1;
            for i = 0.01:0.005:0.5
                gg = @() testCase.example2D(i);
                temps(index) = timeit(gg);
                connecs(index) = size(testCase.example2D(i),1);
                index = index+1;
            end
        end

    end

    methods (Access = private)

        function connec = example2D(testCase, mstep)
            % clc; clear; close all;
%             xdiv = length/mstep +1;
%             ydiv = xdiv
            % Coordinates
            x = 0:  mstep  :1;
            y = 0: mstep/2 :0.25;
            
            %Mesh
            [X,Y,Z] = meshgrid(x,y,0);
            fvc = surf2patch(X,Y,Z,'triangles');
            fvc.vertices(:,3) = []; % 2D
            coords = fvc.vertices;
            connec = fvc.faces;
            m.coord = coords;
            m.connec = connec;
            p.mesh = Mesh(m);
            npnod = size(coords,1);
            
            % Boundary conditions
            dirichNodes = find(coords(:,1) == 0 );
            ndirich = size(dirichNodes,1);
            dirichlet = [dirichNodes,   ones(ndirich,1), zeros(ndirich,1);
                         dirichNodes, 2*ones(ndirich,1), zeros(ndirich,1)];
            neumann   = [npnod, 2, 0.0001];
            
            % FEM parameters
            p.dim = '2D';
            p.type = 'ELASTIC';
            p.scale = 'MACRO';
            p.dirichlet = dirichlet;
            p.pointload = neumann;
            
            % Solution
            fem = NewFEM.create(p);
            fem.solve();
            % fem.plot();
        end
    end

end