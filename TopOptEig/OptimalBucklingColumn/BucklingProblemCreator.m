classdef BucklingProblemCreator < handle

    properties (Access = public)
        designVariable
        cost
        eigenModes
        constraint
    end

    properties (Access = private)
        mesh
        dim
        freeNodes
        inertiaMoment
        youngModulus
    end

    methods
        function obj = BucklingProblemCreator(cParams)
            obj.init(cParams);
            obj.createDesignVariable();
            obj.createCost();
            obj.createEigModes();
            obj.createConstraint();
        end

    end
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh      = cParams.mesh;
            obj.dim       = cParams.dim;
            obj.freeNodes = cParams.freeNodes;
            obj.inertiaMoment = cParams.inertiaMoment;
            obj.youngModulus = cParams.youngModulus;
        end

        function createDesignVariable(obj)
            s.type  = 'AreaColumn';
            s.mesh  = obj.mesh;
            des = DesignVariable.create(s);
            obj.designVariable = des;  
        end

        function createCost(obj)
            sF.type = 'firstEignValue_functional';            
            sC.weights = 1;
            sC.nShapeFuncs = 1;
            sC.designVar = obj.designVariable;         
            sC.shapeFuncSettings{1} = sF;
            obj.cost = Cost(sC);
        end

        function createEigModes(obj)
            s.mesh           = obj.mesh;
            s.inertiaMoment  = obj.inertiaMoment;
            s.youngModulus   = obj.youngModulus;
            s.designVariable = obj.designVariable;
            obj.eigenModes   = EigModes(s);
        end

        function createConstraint(obj)
            sF1.eigModes       = obj.eigenModes;
            sF1.eigNum         = 1;
            sF1.type = 'doubleEig';  

            sF2.eigModes       = obj.eigenModes;  
            sF2.eigNum         = 2;
            sF2.type = 'doubleEig';  

            sF3.type = 'volumeColumn';    
            sF3.mesh = obj.mesh;

            sC.nShapeFuncs = 3;
            sC.designVar = obj.designVariable;   
            sC.dualVariable  = [];
            sC.shapeFuncSettings{1} = sF1;
            sC.shapeFuncSettings{2} = sF2;
            sC.shapeFuncSettings{3} = sF3;
            obj.constraint = Constraint(sC);
        end
    end
end