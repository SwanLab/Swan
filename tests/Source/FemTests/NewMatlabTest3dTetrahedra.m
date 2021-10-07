classdef NewMatlabTest3dTetrahedra < handle
    
    properties (Access = protected)
        testName = 'test3d_tetrahedra';  
        variablesToStore = {'d_u'};
        tol = 1e-6;
        storedVar
        error
        computedVar
        fem
        FileName
    end
       
    methods
        function obj = NewMatlabTest3dTetrahedra()
           obj.computeVariableThroughFemSolver()
           obj.selectComputedVar();
           obj.loadStoredVariable();
        end
    end


    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.fem.variables.d_u;
        end
        
    end
    %% testStoredComputedChecker
    methods (Access = protected) % heredat de testStoredComputedChecker
        function computeError(obj)
            d = numel(obj.variablesToStore);
            err = ones(d,1);
            for ivar = 1:d
                sV = obj.storedVar{ivar};
                cV = obj.computedVar{ivar};
                err(ivar) = norm(sV - cV)/norm(sV);
            end
            obj.error = norm(err);
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

    %% pel suite
    methods (Access = public) % nou
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

