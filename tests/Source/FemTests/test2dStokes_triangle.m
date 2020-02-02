classdef test2dStokes_triangle < testCheckStoredWithFemComputedVariable
    
    
    properties (Access = protected)
        testName = 'test2d_stokes_triangle';  
        variablesToStore = {'variable.u','variable.p'}
    end
       
    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.fem.variables.u; 
            obj.computedVar{2} = obj.fem.variables.p; 
        end  
        
        function computeVariableThroughFemSolver(obj)
            femSolver = FEM.create(obj.testName);
            femSolver.computeVariables;
            obj.fem = femSolver;
        end
       
    end
    

end

