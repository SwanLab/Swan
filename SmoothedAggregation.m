classdef SmoothedAggregation < handle

    properties (Access = private)
        np
        sp
        LHS
        nullSpace
    end

    methods (Access = public)
        function obj = SmoothedAggregation(cParams)
            obj.init(cParams);
        end

        function x = solve(obj,res)

        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.importLibraries();
            obj.setLHS(cParams);
            obj.nullSpace = cParams.nullSpace; % python also
        end

        function importLibraries(obj)
            obj.np = py.importlib.import_module('numpy');
            obj.sp = py.importlib.import_module('scipy.sparse');
        end

        function setLHS(obj,cParams)
            A         = cParams.LHS;
            [i, j, v] = find(A);
            py_i      = obj.np.array(int32(i - 1));
            py_j      = obj.np.array(int32(j - 1));
            py_data   = obj.np.array(v);
            shape_py  = int32(size(A));
            coo       = obj.sp.coo_matrix(py.tuple({py_data, py.tuple({py_i, py_j})}), pyargs('shape', shape_py));
            Acsr      = coo.tocsr();
            obj.LHS   = Acsr;
        end
    end
end