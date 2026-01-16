classdef SmoothedAggregation < handle

    properties (Access = private)
        nLevels
        needUpdate
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
            obj.restartPython();
            obj.init(cParams);
            obj.importLibraries();
            obj.createSolverOptions(cParams);
            obj.setNullSpace(cParams);
        end

        function x = solve(obj,A,res)
            try
                x = obj.compute(A,res);
            catch
                obj.restartPython();
                x = obj.compute(A,res);
            end
        end

        function update(obj)
            obj.needUpdate = true;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.nLevels    = cParams.nLevels;
            obj.needUpdate = true;
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

        function x = compute(obj,A,res)
            if obj.needUpdate
                LHS = obj.convertToPythonSparse(A);
                obj.createAMGSolver(LHS);
                obj.needUpdate = false;
            end
            b    = obj.np.array(double(res));
            x_py = obj.AMGSolver.solve(b, obj.AMGOptions);
            x    = double(x_py)';
        end

        function pyVar = convertToPythonSparse(obj,A)
            [i, j, v] = find(A);
            py_i      = obj.np.array(int32(i - 1));
            py_j      = obj.np.array(int32(j - 1));
            py_data   = obj.np.array(v);

            py_i = py_i.flatten();
            py_j = py_j.flatten();
            py_data = py_data.flatten();

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

    methods (Static, Access = private)
        function restartPython()
            terminate(pyenv);
            pyenv('ExecutionMode','OutOfProcess');
        end
    end
end