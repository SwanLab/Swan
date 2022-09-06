classdef MicroShapeTests < matlab.unittest.TestCase

    properties (TestParameter)
    end

    methods (Test, TestTags = {'MicroShape', 'Hex'})

        function testHexNodeCoords(testCase)
            cases = 'Hex';
            sideLength = [1,1,1];
            theta = [0,60,120];
            divUnit = 3;
            id   = testCase.loadInitialData(cases, sideLength, theta, divUnit);
            t = Tester.create('NodeCoordinatesComputerTester', id);
            res = t.verify();
            testCase.verifyEqual(res, 1);
        end

        function testHexMeshCreator(testCase)
            cases = 'Hex';
            sideLength = [1,1,1];
            theta = [0,60,120];
            divUnit = 3;
            id   = testCase.loadInitialData(cases, sideLength, theta, divUnit);
            t = Tester.create('MeshCreatorTester', id);
            res = t.verify();
            testCase.verifyEqual(res, 1);
            close all;
        end

    end

    methods (Test, TestTags = {'MicroShape', 'Quad'})

        function testQuadNodeCoords(testCase)
            cases = 'Quad';
            sideLength = [1,1];
            theta = [0,90];
            divUnit = 3;
            id   = testCase.loadInitialData(cases, sideLength, theta, divUnit);
            t = Tester.create('NodeCoordinatesComputerTester', id);
            res = t.verify();
            testCase.verifyEqual(res, 1);
        end

        function testQuadMeshCreator(testCase)
            cases = 'Quad';
            sideLength = [1,1];
            theta = [0,90];
            divUnit = 3;
            id   = testCase.loadInitialData(cases, sideLength, theta, divUnit);
            t = Tester.create('MeshCreatorTester', id);
            res = t.verify();
            testCase.verifyEqual(res, 1);
            close all;
        end

    end

    methods (Access = private)

        function initialData = loadInitialData(testCase, cases, sideLength, theta, divUnit)
            nV = load(strcat('nvert',cases,'.mat'));
            bN = load(strcat('boundNodes',cases,'.mat'));
            tN = load(strcat('totalNodes',cases,'.mat'));
            vC = load(strcat('vertCoord',cases,'.mat'));
            b = load(strcat('boundCoord',cases,'.mat'));
            c = load(strcat('totalCoord',cases,'.mat'));
            initialData.c = sideLength;
            initialData.theta = theta;
            initialData.divUnit = divUnit;
            initialData.div = divUnit*sideLength;
            initialData.nvert = nV.nvert;
            initialData.nodes.vert = nV.nvert;
            initialData.nodes.bound = bN.boundNodes;
            initialData.nodes.total = tN.totalNodes;
            initialData.vertCoord = vC.vertCoord;
            initialData.boundCoord = b.boundary;
            initialData.coord = c.coord;
            initialData.filename = 'microshaperesults.m';
        end

    end

end