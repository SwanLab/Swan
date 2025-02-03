classdef BoundaryConditionsTestComputer < handle

    properties (Access = public)
        comparingForceMessage
        comparingDirConditionsMessage
        comparingDirDofsMessage
        comparingNodesMessage
    end

    properties (Access = private)
        
    end

    properties (Access = private)
        forcesFormulaOriginal
        dirConditionsOriginal
        dirDofsOriginal
        nodesConditionsOriginal
        forcesFormulaTest
        dirConditionsTest
        dirDofsTest
        nodesConditionsTest
    end
    
    methods (Access = public)

        function obj = BoundaryConditionsTestComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.compareforcesFormula();
            obj.compareDirConditions();
            obj.compareDirDofs();
            obj.compareNodeConditions();
            result = [obj.comparingForceMessage, newline, obj.comparingDirConditionsMessage, newline, obj.comparingDirDofsMessage, newline, obj.comparingNodesMessage];
            disp(result);
        end

    end


    methods (Access = private)

        function init(obj,cParams)
            obj.saveInput(cParams);
            obj.loadOriginalValues();
        end

        function saveInput(obj,cParams)
            obj.forcesFormulaTest      = cParams.forcesFormula;
            obj.dirConditionsTest      = cParams.dirConditions;
            obj.dirDofsTest            = cParams.dirDofs;
            obj.nodesConditionsTest    = cParams.nodesConditions;
        end

        function loadOriginalValues(obj)
            load("datas.mat","forcesFormula","dirichlet","dir_dofs","nodespresscyl");
            obj.forcesFormulaOriginal   = forcesFormula;
            obj.dirConditionsOriginal   = dirichlet;
            obj.dirDofsOriginal         = dir_dofs;
            obj.nodesConditionsOriginal = nodespresscyl;
        end

        function compareforcesFormula(obj)
           try
                assert(isequal(obj.forcesFormulaOriginal, obj.forcesFormulaTest), ...
                    'The forces formula does not have the expected value.');
                obj.comparingForceMessage = 'The forces formula has the expected value.';
           catch ME
                obj.comparingForceMessage = ME.message;
           end
        end

        function compareDirConditions(obj)
           try
                assert(isequal(obj.dirConditionsOriginal, obj.dirConditionsTest), ...
                    'The dirichlet conditons do not have the expected value.');
                obj.comparingDirConditionsMessage = 'The dirichlet conditons have the expected value.';
           catch ME
                obj.comparingDirConditionsMessage = ME.message;
           end
        end

        function compareDirDofs(obj)
           try
                assert(isequal(obj.dirDofsOriginal, obj.dirDofsTest), ...
                    'The dirichlet dofs do not have the expected value.');
                obj.comparingDirDofsMessage = 'The dirichlet dofs have the expected value.';
           catch ME
                obj.comparingDirDofsMessage = ME.message;
           end
        end

        function compareNodeConditions(obj)
           try
                assert(isequal(obj.nodesConditionsOriginal, obj.nodesConditionsTest), ...
                    'The nodes under the conditions does not match the expected values');
                obj.comparingNodesMessage = 'The nodes under the condition matches the expected values';
           catch ME
                obj.comparingNodesMessage = ME.message;
           end
        end

    end
    
end