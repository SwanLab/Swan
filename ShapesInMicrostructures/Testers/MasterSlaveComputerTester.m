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
            mS = load('masterSlave.mat');
            obj.corrValues(1).Matrix = mS.masterSlave;
        end
        
        function obtainCalculatedData(obj)
            solution = MasterSlaveComputer(obj.data);
            solution.computeMasterSlaveNodes();
            obj.calcValues(1).Matrix = solution.masterSlaveIndex;
        end
    
    end

end