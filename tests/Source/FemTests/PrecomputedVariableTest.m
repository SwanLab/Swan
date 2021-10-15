classdef PrecomputedVariableTest < handle
    
    properties (Access = protected)
        testName;  
        variablesToStore
        storedVar
        computedVar
        computation
        computerType
    end

    methods (Access = public)

        function obj = PrecomputedVariableTest(cParams)
            obj.testName         = cParams.testName;
            obj.variablesToStore = cParams.variablesToStore;
            obj.computerType     = cParams.computerType;
            obj.computeVariable()
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

    end

    methods (Access = private)

        function computeVariable(obj)
            s.testName = obj.testName;
            testComputer = TestComputer.create(obj.computerType, s);
            testComputer.compute();
            obj.computation = testComputer.computation;
        end
        
        function selectComputedVar(obj)
            toStore = obj.variablesToStore;
            comp = obj.computation;
            if isprop(comp, 'variables')
                vars = obj.computation.variables;
            elseif isprop(comp, 'designVariable')
                vars = struct(obj.computation.designVariable);
                valX = toStore{1};
                [vars.(valX)] = vars.value;
                vars = rmfield(vars, 'value');
            end
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
            file2load = obj.testName;
            load(file2load);
            for icell = 1:numel(obj.variablesToStore)
              obj.storedVar{icell} = eval(obj.variablesToStore{icell});
            end
        end

    end

end