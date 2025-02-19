classdef StokesComputer < handle

    properties (Access = public)
        computation
        variables
    end

    properties (Access = private)
        testName
    end

    methods (Access = public)
        function obj = StokesComputer(cParams)
            obj.testName = cParams.testName; % Aquí arribem desde "testComputer = TestComputer.create(obj.computerType, s);", a PrecomputedVariableTest, passant per TestComputer.
        end

        function compute(obj)
            a.fileName = obj.testName;
            s = StokesDataContainer(a); % Anem a "StokesDataContainer" per llegir i emmagatzemar els paràmetres del problema a "s"
            plotejarmalla.plotPunts(s.mesh.coord);
            femSolver = PhysicalProblem.create(s); %Anem a "PhysicalProblem" i d'allà al StokesProblem anem creant els elements per resoldre el problema (com la LHS...)
            femSolver.computeVariables; % Tornem a StokesProblem per resoldre el problema.
            obj.computation = femSolver;
            u = obj.computation.velocityFun.fValues;
            obj.variables.p = obj.computation.pressureFun.fValues;
            obj.variables.u = reshape(u', [numel(u) 1]); % AQUÍ JA TENIM ELS RESULTATS
        end
    end

end

