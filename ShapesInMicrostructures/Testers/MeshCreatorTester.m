classdef MeshCreatorTester < Tester
        
    properties (Access = private)
        data
    end
    
    properties (Access = public) 
        testName
        corrValues
        calcValues
    end    
    
    methods (Access = public)
        
        function obj = MeshCreatorTester(cParams)
            obj.data = cParams;
            obj.testName = 'MeshCreator';
            obj.loadCorrectValues();
            obj.obtainCalculatedData();
        end
        
    end
    
    methods (Access = protected)
        
        function loadCorrectValues(obj)
            switch obj.data.nvert
                case 4
                    c = load('totalCoordQuad.mat');
                    cN = load('connecQuad.mat');
                    mS = load('masterSlaveQuad.mat');
                case 6
                    c = load('totalCoordHex.mat');
                    cN = load('connecHex.mat');
                    mS = load('masterSlaveHex.mat');
            end
            obj.corrValues(1).Matrix = c.coord;
            obj.corrValues(2).Matrix = cN.connec;
            obj.corrValues(3).Matrix = mS.masterSlave;
        end
        
        function obtainCalculatedData(obj)
            solution = MeshCreator(obj.data);
            solution.computeMeshNodes();
            solution.drawMesh();
            solution.plotIndicesOfNodes();
            obj.calcValues(1).Matrix = solution.coord;
            obj.calcValues(2).Matrix = solution.connec;
            obj.calcValues(3).Matrix = solution.masterSlaveIndex;
        end
    
    end
      
end