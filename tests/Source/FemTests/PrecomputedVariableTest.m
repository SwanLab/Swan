classdef PrecomputedVariableTest < handle
    
    properties (Access = protected)
        testName;  
        variablesToStore
        storedVar
        computedVar
        fem
        FileName
        TestComputationHandler
        TestVariableHandler
    end
       
    methods
        function obj = PrecomputedVariableTest(cParams)
            obj.testName         = cParams.testName;
            obj.variablesToStore = cParams.variablesToStore;
            obj.computeVariableThroughFemSolver()
            obj.selectComputedVar();
            obj.loadStoredVariable();
        end
    end


   
    methods (Access = protected)

        function computeVariableThroughFemSolver(obj)
            obj.fem = FEM.create(obj.testName);
            obj.createMaterialProperties();
            obj.fem.computeVariables();
        end
        
        function selectComputedVar(obj)
            vars = obj.fem.variables;
            toStore = obj.variablesToStore;
            fnms = fieldnames(vars);
            totalComputationVariables = numel(fnms);
            totalStoredVariables = size(toStore,2);
            count = 1;
            for i = 1:totalComputationVariables
                calcVar = fnms(i);
                storVar = toStore(count);
                if strcmp(calcVar, storVar)
                    obj.computedVar{count} = vars.(fnms{i});
                    if count == totalStoredVariables
                        return;
                    else
                        count = count+1;
                    end
                end

            end
        end
        
    end

    %% testLoadStoredVariable
    
    methods (Access = private)
        
        function loadStoredVariable(obj) % heredat de testLoadStoredVariable
            file2load = obj.testName;
            load(file2load);
            for icell = 1:numel(obj.variablesToStore)
              obj.storedVar{icell} = eval(obj.variablesToStore{icell});
            end
        end

        function createMaterialProperties(obj) % heredat de testFemComputation
            q = Quadrature.set(obj.fem.mesh.type);
            q.computeQuadrature('LINEAR');
            I = ones(obj.fem.mesh.nelem,q.ngaus);
            p.kappa = .9107*I;
            p.mu    = .3446*I;
            obj.fem.setMatProps(p)     
        end

    end

    methods (Access = public)

        function error =computeError(obj)
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

end

