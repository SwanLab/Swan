classdef MasterSlaveComputerTester < Tester
            
    properties (Access = private)
        data
    end
    
    properties (Access = protected) 
        testName
        corrValues
        calcValues
    end    
    
    methods (Access = public)
        
        function obj = MasterSlaveComputerTester(cParams)
            obj.data = cParams;
            obj.testName = 'MasterSlaveComputerTester';
            obj.loadCorrectValues();
            obj.obtainCalculatedData();
            obj.verify();
        end
        
    end
    
    methods (Access = protected)
        
        function loadCorrectValues(obj)
            switch obj.data.nvert
                case 4
                    mS = load('masterSlaveQuad.mat');
                case 6
                    mS = load('masterSlaveHex.mat');
            end
            obj.corrValues(1).Matrix = mS.masterSlave;
        end
        
        function obtainCalculatedData(obj)
            solution = MasterSlaveComputer(obj.data);
            solution.computeMasterSlaveNodes();
            obj.calcValues(1).Matrix = solution.masterSlaveIndex;
        end
    
    end

end