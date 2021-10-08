classdef NewMatlabTest < handle
    
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
        function obj = NewMatlabTest(cParams)
            obj.testName         = cParams.testName;
            obj.variablesToStore = cParams.variablesToStore;
%             solverType = cParams.solver;
%             obj.TestComputationHandler = NewTestComputationHandler.create(solverType);

            obj.computeVariableThroughFemSolver()
            obj.selectComputedVar();
            obj.loadStoredVariable();
        end
    end


    methods (Access = protected)
        
        function selectComputedVar(obj)
            vars = obj.fem.variables;
            toStore = obj.variablesToStore;
            fnms = fieldnames(vars);
            count = 1;
            for i = 1:numel(vars)
                if strcmp(fnms(i),toStore{i})
                    obj.computedVar{count} = vars.(fnms{i});
                    count = count + 1;
                end
            end
        end
        
    end

    %% testLoadStoredVariable
    
    methods (Access = private) % heredat de testLoadStoredVariable
        function loadStoredVariable(obj)
            file2load = obj.testName;
            load(file2load);
            for icell = 1:numel(obj.variablesToStore)
              obj.storedVar{icell} = eval(obj.variablesToStore{icell});
            end
        end
    end

    %% testFemComputation
    
    methods (Access = protected) % heredat de testFemComputation
        function computeVariableThroughFemSolver(obj)
            obj.fem = FEM.create(obj.testName);
            obj.createMaterialProperties();
            obj.fem.computeVariables();
        end        
    end
    
    methods (Access = private) % heredat de testFemComputation
        
        function createMaterialProperties(obj)
            q = Quadrature.set(obj.fem.mesh.type);
            q.computeQuadrature('LINEAR');
            I = ones(obj.fem.mesh.nelem,q.ngaus);
            p.kappa = .9107*I;
            p.mu    = .3446*I;
            obj.fem.setMatProps(p)     
        end

    end

    methods (Access = public)

        function error =computeErrorForTest(obj)
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

