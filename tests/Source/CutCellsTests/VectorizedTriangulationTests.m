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
                             'OneTetrahedronIsoCoordOrderedConnecCase7';
                            'OneTetrahedronRandCoordOrderedConnecCase7';
%                            'OneTetrahedronIsoCoordOrderedConnecCase11';
%                           'OneTetrahedronRandCoordOrderedConnecCase11';
%                            'OneTetrahedronIsoCoordOrderedConnecCase13';
%                           'OneTetrahedronRandCoordOrderedConnecCase13';
%                            'OneTetrahedronIsoCoordOrderedConnecCase8';
%                           'OneTetrahedronRandCoordOrderedConnecCase8';
%                             'OneTetrahedronIsoCoordOrderedConnecCase4';
%                           'OneTetrahedronRandCoordOrderedConnecCase4';
%                             'OneTetrahedronIsoCoordOrderedConnecCase2';
%                           'OneTetrahedronRandCoordOrderedConnecCase2';
%                             'OneTetrahedronIsoCoordOrderedConnecCase14';
%                           'OneTetrahedronRandCoordOrderedConnecCase14';
%                             'OneTetrahedronIsoCoordOrderedConnecCase1';
%                           'OneTetrahedronRandCoordOrderedConnecCase1';
%                            'OneTetrahedronAllRand3Vs1';
%                              'TwoTetrahedronIsoCoordOrderedConnecCase3VsOnePoint'
%                              'OneTetrahedronIsoCoordOrderedConnecCase6';
%                              'OneTetrahedronIsoCoordOrderedConnecCase9';
%                               'OneTetrahedronIsoCoordOrderedConnecCase12';
%                               'OneTetrahedronIsoCoordOrderedConnecCase3';
%                               'OneTetrahedronIsoCoordOrderedConnecCase10';
%                               'OneTetrahedronIsoCoordOrderedConnecCase5';
%                               'OneTetrahedronAllRand2Vs2';
%                              'TwoTetrahedronIsoCoordOrderedConnecCase2Vs2Point';
%                               'TwoTetrahedronIsoCoordOrderedConnecCaseBothCases';
%               'TwoTetrahedronIsoCoordRandConnecCaseBothCases';  
%               'FindingBug'
            };
        end
    end
end
