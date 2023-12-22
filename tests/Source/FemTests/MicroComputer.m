classdef MicroComputer < handle

    properties (Access = public)
        computation
        variables
    end

    properties (Access = private)
        testName
        builderType
        solType
        solMode
    end

    methods (Access = public)
        function obj = MicroComputer(cParams)
            obj.testName     = cParams.testName;
%             obj.builderType  = cParams.builderType;
            obj.solType = cParams.solType;
            obj.solMode = cParams.solMode;
        end

        function compute(obj)
            a.fileName = obj.testName;
            s = FemDataContainer(a);

            % s = ObjectSaver.saveObj(s);
%             s.builderType = obj.builderType;
            s.solType = obj.solType;
            s.solMode = obj.solMode;
            switch s.solMode
                case 'DISP'
                    femSolver = ElasticProblemDisp(s);
                    femSolver.computeStressHomog();
                    obj.computation = femSolver;
                case 'FLUC'
                    femSolver = ElasticProblemFluc(s);
                    switch s.solType
                        case 'REDUCED'
                            femSolver.computeChomog();
                            obj.computation = femSolver;
                        case 'MONOLITIC'
                            femSolver.computeLagrangeMultSum();
                            obj.computation = femSolver;
                    end

            end
%             femSolver = ElasticProblemMicro(s);
%             femSolver.computeChomog();
%             obj.computation = femSolver;

%             femSolver = ElasticProblemMicro(s);
%             femSolver = NewElasticProblemMicro(s);
            femSolver = ElasticProblemMicro(s);
            femSolver.solve();
            obj.computation = femSolver;
            obj.variables.Chomog = femSolver.Chomog;
        end

    end

end
