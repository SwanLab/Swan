classdef TestOneTetrahedronAllRand <  handle
    
    methods (Static)

        function test = create(test_type)
            s.coord = rand(4,3);
            s.connec = randperm(4);
            s.boundaryConnec = [1 2 3;1 2 4;1 3 4;2 3 4];
            s.connecBcutMesh = [5 6 7]; % afegit meu
            b = 10;
            a = 0;
            ls = rand(4,1);
            ls = a + (b-a)*ls;
            switch test_type
                case {'2Vs2'}
                    position = randperm(4,2);
                    ls(position) = -ls(position);
                    s.levelSet = ls;
                    test = VectorizedTriangulationTest(s);
                case {'3Vs1'}
                    difPos = randperm(4,1);
                    isPositive = randperm(2,1);
                    position = false(4,1);
                    if isPositive == 1
                        position(difPos) = true;
                    else
                        position(setdiff(1:4,difPos)) = true;
                    end
                    ls(position) = -ls(position);
                    s.levelSet = ls;
                    test = VectorizedTriangulationTest(s);
                otherwise
                    error('Invalid Tetrahedron Test Type.')
            end
        end

    end

    methods (Static)

        function s = createCoordAndConnec()
            s.coord = rand(4,3);
            s.connec = randperm(4);
            s.boundaryConnec = [1 2 3;1 2 4;1 3 4;2 3 4];
        end

        function createLevelSet2vs2(s)
            b = 10;
            a = 0;
            ls = rand(4,1);
            ls = a + (b-a)*ls;
            position = randperm(4,2);
            ls(position) = -ls(position);
            s.levelSet = ls;
        end

        function createLevelSet3vs1(s)
            b = 10;
            a = 0;
            ls = rand(4,1);
            ls = a + (b-a)*ls;
            difPos = randperm(4,1);
            isPositive = randperm(2,1);
            position = false(4,1);
            if isPositive == 1
                position(difPos) = true;
            else
                position(setdiff(1:4,difPos)) = true;
            end
            ls(position) = -ls(position);
            s.levelSet = ls;
        end

    end

end