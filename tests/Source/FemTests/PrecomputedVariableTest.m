classdef PrecomputedVariableTest < handle

    properties (Access = protected)
        testName  
        variablesToStore
        storedVar
        computedVar
        computation
        computerType
        testResultsName
        computer
    end

%% Mètodes pblics
    methods (Access = public)

        function obj = PrecomputedVariableTest(cParams) %cParams és la estructura s
            obj.init(cParams); % Anem a la funció "init" ⬇
            obj.computeVariable(cParams); % Anem a la resta de funcions . . .
            obj.selectComputedVar();
            obj.loadStoredVariable();
        end

        function error = computeError(obj)
            d = numel(obj.variablesToStore);
            err = ones(d,1);
            for ivar = 1:d
                sV = obj.storedVar{ivar};
                cV = obj.computedVar{ivar};
                err(ivar) = norm(sV - cV)/norm(sV);
            end
            error = norm(err);
        end

        function overwriteResults(obj)
            % NOTE: change the variable name "x" to the one in the .mat
            % file so that it works
            variableNames = who('-file', obj.testResultsName);
            d = numel(variableNames);
            for ivar = 1:d
                varName = variableNames{ivar};
                eval([varName ' = obj.computedVar{ivar}' ])
                save(obj.testResultsName, varName, '-append')
            end
        end

        function comp = getComputation(obj)
            comp = obj.computation;
        end
    end

%% Mètodes privats

    methods (Access = private)

        function init(obj, cParams) % Asignem els paràmetres:
            obj.testName         = cParams.testName; 
            obj.variablesToStore = cParams.variablesToStore;
            obj.computerType     = cParams.computerType;
            obj.testResultsName  = cParams.testName;
            if isfield(cParams, 'testResultsName')
                obj.testResultsName = cParams.testResultsName;
            end
        end

        function computeVariable(obj, s)
            s.testName = obj.testName;
            testComputer = TestComputer.create(obj.computerType, s); % Anem a "TestComputer" que ens portarà a "StokesComputer"
            testComputer.compute(); % Ara li demanem a la classe "StokesComputer" que utilizi el seu mètode "compute"
            obj.computation = testComputer.computation;
            obj.computer = testComputer;
        end
        
        function selectComputedVar(obj)
            toStore = obj.variablesToStore;
            vars = obj.computer.variables;
            fnms = fieldnames(vars);
            totalComputationVariables = numel(fnms);
            totalStoredVariables = size(toStore,2);
            for i = 1:totalComputationVariables
                calcVar = fnms(i);
                for j = 1:totalStoredVariables
                    storVar = toStore(j);
                    if strcmp(calcVar, storVar)
                        obj.computedVar{j} = vars.(fnms{i});
                        break
                    end
                end
            end
        end

        function loadStoredVariable(obj)
            file2load = obj.testResultsName;
            load(file2load);
            for icell = 1:numel(obj.variablesToStore)
              obj.storedVar{icell} = eval(obj.variablesToStore{icell});
            end
        end

    end

end