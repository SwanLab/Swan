classdef VectorizedTriangulationTests < testRunner
    properties (Access = protected)
        FieldOfStudy = 'Vectorized Triangulation'
        tests
    end
    
    methods (Access = public)
        function  obj = VectorizedTriangulationTests()
            obj@testRunner();
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {... 
%                             'OneTetrahedronRandCoordOrderedConnecCase7'; %funciona
%                             'OneTetrahedronRandCoordOrderedConnecCase11';
%                             'OneTetrahedronRandCoordOrderedConnecCase13';
%                             'OneTetrahedronRandCoordOrderedConnecCase8';
%                             'OneTetrahedronRandCoordOrderedConnecCase4';
%                             'OneTetrahedronRandCoordOrderedConnecCase2';
%                             'OneTetrahedronRandCoordOrderedConnecCase14';
%                             'OneTetrahedronRandCoordOrderedConnecCase1';
%                               'OneTetrahedronIsoCoordOrderedConnecCase1';
%                               'OneTetrahedronIsoCoordOrderedConnecCase2';
%                               'OneTetrahedronIsoCoordOrderedConnecCase3';
%                               'OneTetrahedronIsoCoordOrderedConnecCase4';
%                               'OneTetrahedronIsoCoordOrderedConnecCase5';
%                               'OneTetrahedronIsoCoordOrderedConnecCase6';
%                               'OneTetrahedronIsoCoordOrderedConnecCase7'; %funciona
%                               'OneTetrahedronIsoCoordOrderedConnecCase8';
%                               'OneTetrahedronIsoCoordOrderedConnecCase9';
%                               'OneTetrahedronIsoCoordOrderedConnecCase10';
%                               'OneTetrahedronIsoCoordOrderedConnecCase11';
%                               'OneTetrahedronIsoCoordOrderedConnecCase12';
%                               'OneTetrahedronIsoCoordOrderedConnecCase13';
%                               'OneTetrahedronIsoCoordOrderedConnecCase14';
%               'FindingBug'
%               'OneTetrahedronAllRand3Vs1';
%               'OneTetrahedronAllRand2Vs2';
%               'TwoTetrahedronIsoCoordOrderedConnecCase3VsOnePoint'
%               'TwoTetrahedronIsoCoordOrderedConnecCase2Vs2Point';
%               'TwoTetrahedronIsoCoordOrderedConnecCaseBothCases';
%               'TwoTetrahedronIsoCoordRandConnecCaseBothCases';  
            };
        end
    end
end
