classdef NewVectorizedTriangulationTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        levelSet1 = struct('positive', 7.8496, 'negative', -7.8496);
        levelSet2 = struct('positive', 9.7731, 'negative', -9.7731);
        levelSet3 = struct('positive', 8.3404, 'negative', -8.3404);
        lvlSet4 = struct('positive', 8.3622, 'negative', -8.3622);
        allRandTests = {'2Vs2', '3Vs1'}
    end

    methods (Test, ...
            TestTags = {'VectorizedTriangulation', 'Classic', 'Tetrahedron', 'IsoCoord', 'OrderedConnec'})
        function testIsoCoord(testCase, levelSet1, levelSet2, levelSet3, lvlSet4)
            s.coord    = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
            s.connec   = [1 2 3 4];
            s.levelSet = [levelSet1; levelSet2; levelSet3; lvlSet4];
            s.boundaryConnec = [1 2 3;1 2 4;1 3 4;2 3 4];
            s.connecBcutMesh = [5 6 7];
            test = NewVectorizedTriangulationTest(s);
            err = test.computeError();
            tol = 1e-14;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, ...
            TestTags = {'VectorizedTriangulation', 'Classic', 'Tetrahedron', 'IsoCoord', 'OrderedConnec'})
        function testRandCoord(testCase, levelSet1, levelSet2, levelSet3, lvlSet4)
            s.coord    = rand(4,3);
            s.connec   = [1 2 3 4];
            s.levelSet = [levelSet1; levelSet2; levelSet3; lvlSet4];
            s.boundaryConnec = [1 2 3;1 2 4;1 3 4;2 3 4];
            s.connecBcutMesh = [5 6 7];
            test = NewVectorizedTriangulationTest(s);
            err = test.computeError();
            tol = 1e-14;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, ...
            TestTags = {'VectorizedTriangulation', 'Classic', 'Tetrahedron', 'IsoCoord', 'OrderedConnec'})
        function testFindingBug(testCase)
            s.coord    = [0    0    0;
                          1    0    0;
                          0    1    0;
                          0.1  0.1  0.7];
            s.connec   = [2 1 3 4];
            s.levelSet = [-3.9507; 15.2486; -10.0641; 5.1921];
            s.boundaryConnec = [2 4 1;2 4 3;2 1 3;4 1 3];
            s.connecBcutMesh = [5 6 7]; % afegit meu
            test = NewVectorizedTriangulationTest(s);
            err = test.computeError();
            tol = 1e-14;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, ...
            TestTags = {'VectorizedTriangulation', 'Classic', 'AllRand'})

        function testOneTetrahedronAllRand (testCase,allRandTests)
            test = NewTestOneTetrahedronAllRand.create(allRandTests);
            err = test.computeError();
            tol = 1e-14;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, ...
            TestTags = {'VectorizedTriangulation', 'Nou', 'TwoTetrahedron'})
        
        function testTwoTetrahedron2Vs2Point (testCase)
            s.coord = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 1];
            s.connec = [1 2 3 4;
                        2 3 4 5];
            s.boundaryConnec = [1 2 3;1 2 4;1 3 4;5 2 3;5 2 4;5 3 4];
            s.levelSet =  [7.8496; 5.7731; -8.3404; -8.3622; 10.3622];
            s.connecBcutMesh = [5 6 7]; % collita propia
            test = NewVectorizedTriangulationTest(s);
            err = test.computeError();
            tol = 1e-14;
            testCase.verifyLessThanOrEqual(err, tol)
        end
        
        function testTwoTetrahedron3VsOnePoint (testCase)
            s.coord = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 1];
            s.connec = [1 2 3 4;
                        2 3 4 5];
            s.boundaryConnec = [1 2 3; 1 2 4; 1 3 4; 5 2 3; 5 2 4; 5 3 4];
            s.levelSet = [-7.8496; 5.7731; -8.3404; -8.3622; -10.3622];
            s.connecBcutMesh = [5 6 7]; % collita propia
            test = NewVectorizedTriangulationTest(s);
            err = test.computeError();
            tol = 1e-14;
            testCase.verifyLessThanOrEqual(err, tol)
        end
        
        function testTwoTetrahedronOrderedBothCases (testCase)
            s.coord = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 1];
            s.connec = [1 2 3 4;
                        2 3 4 5];
            s.boundaryConnec = [1 2 3;1 2 4;1 3 4;5 2 3;5 2 4;5 3 4];
            s.connecBcutMesh = [6 5 7]; % afegit meu
            s.levelSet = [-7.8496; 5.7731; -8.3404; -8.3622; 10.3622];
            test = NewVectorizedTriangulationTest(s);
            err = test.computeError();
            tol = 1e-14;
            testCase.verifyLessThanOrEqual(err, tol)
        end
        
        function testTwoTetrahedronRandBothCases (testCase)
            s.coord = rand(50,3); 
            t = delaunayTriangulation(s.coord);
            s.connec = t.ConnectivityList;
            s.boundaryConnec = boundary(s.coord);
            ls = 10 * rand(size(s.coord,1),1);
            number = randperm(size(s.coord,1),1);
            position = randperm(size(s.coord,1),number);      
            ls(position) = -ls(position);
            s.levelSet = ls;
            s.connecBcutMesh = [6 5 7]; % afegit meu
            test = NewVectorizedTriangulationTest(s);
            err = test.computeError();
            tol = 1e-14;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

end