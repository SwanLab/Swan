classdef VectorizedTriangulationTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        levelSet1 = struct('pos', 7.8496, 'neg', -7.8496);
        levelSet2 = struct('pos', 9.7731, 'neg', -9.7731);
        levelSet3 = struct('pos', 8.3404, 'neg', -8.3404);
        levelSet4 = struct('pos', 8.3622, 'neg', -8.3622);
        allRandTests = {'2Vs2', '3Vs1'}
    end

    methods (Test, ...
            TestTags = {'VectorizedTriangulation', ...
            'Tetrahedron', 'IsoCoord', 'OrderedConnec'})
        function testIsoCoord(testCase, levelSet1, levelSet2, levelSet3, levelSet4)
            sign1 = sign(levelSet1);
            sign2 = sign(levelSet2);
            sign3 = sign(levelSet3);
            sign4 = sign(levelSet4);
            paramCombination = [sign1, sign2, sign3, sign4];
            invalid = [
                    +1 +1 +1 +1;
                    +1 +1 -1 -1;
                    +1 -1 +1 -1;
                    +1 -1 -1 +1; %
                    -1 +1 +1 -1;
                    -1 +1 -1 +1;
                    -1 -1 +1 +1;
                    -1 -1 -1 -1;
                ];
            isItValidTest = ~ismember(paramCombination,invalid,'rows');

            if (isItValidTest)
                s.coord    = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
                s.connec   = [1 2 3 4];
                s.levelSet = [levelSet1; levelSet2; levelSet3; levelSet4];
                s.boundaryConnec = [1 2 3;1 2 4;1 3 4;2 3 4];
                s.connecBcutMesh = [1 2 3];
                test = VectorizedTriangulationTest(s);
                err = test.computeError();
                tol = 1e-14;
                testCase.verifyLessThanOrEqual(err, tol)
            else
                passed = true;
                verifyTrue(testCase, passed)
            end
        end

    end

    methods (Test, ...
            TestTags = {'VectorizedTriangulation', ...
            'Tetrahedron', 'RandCoord', 'OrderedConnec'})
        function testRandCoord(testCase, levelSet1, levelSet2, levelSet3, levelSet4)
            sign1 = sign(levelSet1);
            sign2 = sign(levelSet2);
            sign3 = sign(levelSet3);
            sign4 = sign(levelSet4);
            paramCombination = [sign1, sign2, sign3, sign4];
            valid = [
                    -1 +1 +1 +1;
                    +1 -1 +1 +1;
                    +1 +1 -1 +1;
                    -1 -1 -1 +1;
                    +1 +1 +1 -1;
                    -1 -1 +1 -1;
                    -1 +1 -1 -1;
                    +1 -1 -1 -1;
                ];
            isItValidTest = ismember(paramCombination,valid,'rows');

            if (isItValidTest)
                s.coord    = rand(4,3);
                s.connec   = [1 2 3 4];
                s.levelSet = [levelSet1; levelSet2; levelSet3; levelSet4];
                s.boundaryConnec = [1 2 3;1 2 4;1 3 4;2 3 4];
                s.connecBcutMesh = [1 2 3];
                test = VectorizedTriangulationTest(s);
                err = test.computeError();
                tol = 1e-14;
                testCase.verifyLessThanOrEqual(err, tol)
            else
                passed = true;
                verifyTrue(testCase, passed)
            end
        end

    end

    methods (Test, ...
            TestTags = {'VectorizedTriangulation', ...
            'Tetrahedron', 'FindingBug', 'OrderedConnec', 'ShallPass'})
        function testFindingBug(testCase)
            s.coord    = [0    0    0;
                          1    0    0;
                          0    1    0;
                          0.1  0.1  0.7];
            s.connec   = [2 1 3 4];
            s.levelSet = [-3.9507; 15.2486; -10.0641; 5.1921];
            s.boundaryConnec = [2 4 1;2 4 3;2 1 3;4 1 3];
            s.connecBcutMesh = [1 2 3]; % Irrelevant but needed
            test = VectorizedTriangulationTest(s);
            passed = true;
            verifyTrue(testCase, passed)
        end

    end

    methods (Test, ...
            TestTags = {'VectorizedTriangulation', 'ShallPass', 'AllRand'})

        function testOneTetrahedronAllRand (testCase,allRandTests)
            test = TestOneTetrahedronAllRand.create(allRandTests);
            passed = true;
            verifyTrue(testCase, passed)
        end

    end

    methods (Test, ...
            TestTags = {'VectorizedTriangulation', 'ShallPass', 'TwoTetrahedron'})
        
        function testTwoTetrahedron2Vs2Point (testCase)
            s.coord = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 1];
            s.connec = [1 2 3 4;
                        2 3 4 5];
            s.boundaryConnec = [1 2 3;1 2 4;1 3 4;5 2 3;5 2 4;5 3 4];
            s.levelSet =  [7.8496; 5.7731; -8.3404; -8.3622; 10.3622];
            s.connecBcutMesh = [1 2 3]; % Irrelevant but needed
            test = VectorizedTriangulationTest(s);
            passed = true;
            verifyTrue(testCase, passed)
        end
        
        function testTwoTetrahedron3VsOnePoint (testCase)
            s.coord = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 1];
            s.connec = [1 2 3 4;
                        2 3 4 5];
            s.boundaryConnec = [1 2 3; 1 2 4; 1 3 4; 5 2 3; 5 2 4; 5 3 4];
            s.levelSet = [-7.8496; 5.7731; -8.3404; -8.3622; -10.3622];
            s.connecBcutMesh = [1 2 3]; % Irrelevant but needed
            test = VectorizedTriangulationTest(s);
            passed = true;
            verifyTrue(testCase, passed)
        end
        
        function testTwoTetrahedronOrderedBothCases (testCase)
            s.coord = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 1];
            s.connec = [1 2 3 4;
                        2 3 4 5];
            s.boundaryConnec = [1 2 3;1 2 4;1 3 4;5 2 3;5 2 4;5 3 4];
            s.connecBcutMesh = [1 2 3]; % Irrelevant but needed
            s.levelSet = [-7.8496; 5.7731; -8.3404; -8.3622; 10.3622];
            test = VectorizedTriangulationTest(s);
            passed = true;
            verifyTrue(testCase, passed)
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
            s.connecBcutMesh = [1 2 3]; % Irrelevant but needed
            test = VectorizedTriangulationTest(s);
            passed = true;
            verifyTrue(testCase, passed)
        end

    end

end