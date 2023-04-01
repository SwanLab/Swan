classdef CostComputer < handle
    properties (Access = public)
        cost
        derivedCost
    end
    properties (Access = private)
        mesh
        structure
        designFields
        displacement
    end 
    methods (Access = public)
        function obj = CostComputer(cParams)
            obj.inputData(cParams)            
        end
        function computeCost(obj)
            s.mesh = obj.mesh;
            s.structure =obj.structure;
            s.projectedField = obj.designFields.projectedField;
            B = FEMcomputer(s);
            B.compute();
            obj.displacement = B.displacement;
            force = B.force;
            obj.cost  = abs(force'*obj.displacement);
        end
        function deriveCost(obj)
            s.mesh = obj.mesh;
            s.displacement =obj.displacement;
            s.structure =obj.structure;
            s.designFields = obj.designFields;
            s.cost =obj.cost ;    
            B = CostFieldDerivator(s);
            B.compute();
            obj.derivedCost = B.derivedCost;
        end 
    end
    methods (Access = private)
        function inputData(obj,cParams) 
            obj.mesh = cParams.mesh;
            obj.structure =cParams.structure;
            obj.designFields = cParams.designFields;
        end      
    end 
end