classdef SmoothedAggregation < handle

    properties (Access = private)
        nLevels
        nSolve
        np
        sp
        AMGOptions
        nullSpace
    end

    properties (Access = private)
        AMGSolver
    end

    methods (Access = public)
        function obj = SmoothedAggregation(cParams)
            obj.init(cParams);
            obj.importLibraries();
            obj.createSolverOptions(cParams);
            obj.setNullSpace(cParams);
        end

        function x = solve(obj,A,res)
            if obj.nSolve/50 == round(obj.nSolve/50)
                LHS = obj.convertToPythonSparse(A);
                obj.createAMGSolver(LHS);
            end
            b    = obj.np.array(double(res));
            x_py = obj.AMGSolver.solve(b, obj.AMGOptions);
            x    = double(x_py)';
            obj.nSolve = obj.nSolve + 1;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.nLevels = cParams.nLevels;
            obj.nSolve  = 0;
        end

        function importLibraries(obj)
            obj.np = py.importlib.import_module('numpy');
            obj.sp = py.importlib.import_module('scipy.sparse');
        end

        function createSolverOptions(obj,cParams)
            tol     = cParams.tol;
            maxIter = cParams.maxIter;
            obj.AMGOptions = pyargs("tol", tol, "maxiter", maxIter);
        end

        function setNullSpace(obj,cParams)
            BNull         = cParams.nullSpace;
            obj.nullSpace = obj.np.array(full(BNull));
        end

        function pyVar = convertToPythonSparse(obj,A)
            [i, j, v] = find(A);
            py_i      = obj.np.array(int32(i - 1));
            py_j      = obj.np.array(int32(j - 1));
            py_data   = obj.np.array(v);
            shape_py  = int32(size(A));
            coo       = obj.sp.coo_matrix(py.tuple({py_data, py.tuple({py_i, py_j})}), pyargs('shape', shape_py));
            pyVar     = coo.tocsr();
        end

        function createAMGSolver(obj,LHS)
            pyArgs        = pyargs('B', obj.nullSpace,'max_levels', int32(obj.nLevels));
            s             = {LHS};
            s{end+1}      = pyArgs;
            obj.AMGSolver = py.pyamg.smoothed_aggregation_solver(s{:});
        end
    end
end