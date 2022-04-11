classdef MasterSlaveNodesTester < Tester
    
    properties (Access = private)
        data
    end
    
    properties (Access = protected)
        testName
        corrValues
        calcValues
    end
    
    methods (Access = public)
        
        function obj = MasterSlaveNodesTester(cParams)
            obj.data = cParams;
            obj.testName = 'MasterSlaveNodesComputer';
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
            solution = computeMasterSlaveNodes(obj.data.vert,obj.data.bound,obj.data.nsides,obj.data.div,obj.data.dim);
            obj.calcValues(1).Matrix = solution;
        end
    
    end
    
end