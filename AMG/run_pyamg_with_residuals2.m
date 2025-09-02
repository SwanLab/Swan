function [x, residuals] = run_pyamg_with_residuals2(A, b, varargin)
%RUN_PYAMG_WITH_RESIDUALS Solve Ax = b using PyAMG with residual tracking.
%
% Input:
%   A - MATLAB sparse matrix
%   b - RHS vector (column vector)
%   varargin - Optional:
%              'B_nullspace', Near-nullspace modes (NumPy array or MATLAB matrix)
%
% Output:
%   x - Solution vector
%   residuals - Residual norms at each iteration

    % --- Parse Optional Inputs ---
    p = inputParser;
    addParameter(p, 'B_nullspace', [], @(x) isempty(x) || isnumeric(x)); % Allow empty or numeric
    parse(p, varargin{:});

    B_nullspace_matlab = p.Results.B_nullspace;
    % --- End Optional Input Parsing ---

    % --- Hardcoded Solver Parameters ---
    % These values are now fixed within the function
    solver_tol = 1e-8;    % Default tolerance for the solver
    solver_maxiter = 1000; % Default maximum iterations for the solver
    % --- End Hardcoded Solver Parameters ---

    % --- Input Validation (MATLAB side) ---
    % Ensure A is a square matrix for typical Ax=b problems
    if size(A, 1) ~= size(A, 2)
        error('Input matrix A must be square (number of rows must equal number of columns).');
    end
    % Ensure the number of rows in A matches the number of elements in b
    if size(A, 1) ~= numel(b)
        error('Dimension mismatch: The number of rows in A must match the number of elements in b.');
    end
    % If B_nullspace is provided, ensure its row dimension matches A's
    if ~isempty(B_nullspace_matlab) && size(B_nullspace_matlab, 1) ~= size(A, 1)
        error('Dimension mismatch: The number of rows in B_nullspace must match the number of rows in A.');
    end
    % --- End Input Validation ---

    % Ensure Python is initialized
    if ~strcmp(pyenv().Status, "Loaded")
        error("Python must be initialized with pyenv first.");
    end

    % Define the Python file name
    pyfile = 'ResidualTracker.py';

    % --- Explicitly delete the file if it exists, to ensure a fresh copy ---
    if isfile(pyfile)
        try
            delete(pyfile);
            fprintf('Deleted existing %s.\n', pyfile);
        catch ME
            warning('Could not delete existing %s: %s', pyfile, ME.message);
        end
    end
    % --- End delete file ---

    % Write Python callback class to file
    fid = fopen(pyfile, 'w');
    if fid == -1
        error('Could not open %s for writing. Check file permissions or directory.', pyfile);
    end
    fprintf(fid, '%s\n', ...
        "import numpy as np", ... % Import numpy for norm calculation
          "class ResidualTracker:", ...
          "    def __init__(self, A, b):", ...
          "        self.residuals = []", ...
          "        self.A = A", ...
          "        self.b = b", ...
          "", ...
          "    def __call__(self, xk):", ...
          "        r = self.b - self.A @ xk", ...
          "        norm_r = np.linalg.norm(r)", ...
          "        self.residuals.append(float(norm_r))");
    fclose(fid);

    % Add current folder to Python path if not already present
    % This ensures MATLAB can find 'ResidualTracker.py'
    current_folder = fileparts(mfilename('fullpath'));
    if ~any(cellfun(@(p) strcmp(p, current_folder), cell(py.sys.path)))
        insert(py.sys.path, int32(0), current_folder);
    end

    % Import required Python modules
    % Make sure pyamg is installed in your Python environment (e.g., pip install pyamg)
    py.importlib.import_module('pyamg');
    ResidualTracker = py.importlib.import_module('ResidualTracker');
    np = py.importlib.import_module('numpy');
    sp = py.importlib.import_module('scipy.sparse');

    % Convert MATLAB sparse matrix A to Python scipy.sparse.csr_matrix
    [i, j, v] = find(A);

    % Convert to zero-based indexing for Python and NumPy arrays
    py_i = np.array(int32(i - 1));
    py_j = np.array(int32(j - 1));
    py_data = np.array(v);

    % Get the shape of the matrix
    shape_py = int32(size(A));

    % Correctly construct the coo_matrix.
    % The data, row, and col are passed as a tuple: (data, (row, col))
    % The shape is passed as a keyword argument.
    coo = sp.coo_matrix(py.tuple({py_data, py.tuple({py_i, py_j})}), pyargs('shape', shape_py));
    A_csr = coo.tocsr();

    % Convert RHS b to numpy array (ensure it's a double and column vector)
    b_np = np.array(double(b(:)));

    % Create residual tracker instance
%     tracker = ResidualTracker.ResidualTracker();
    tracker = ResidualTracker.ResidualTracker(A_csr, b_np);
    residuals = [];
    % Define standard solver options for the solve method
%     solve_method_options = pyargs("callback", tracker, "tol", solver_tol, "maxiter", solver_maxiter);
    solve_method_options = pyargs("callback", tracker, "tol", solver_tol, "maxiter", solver_maxiter, "accel", "cg");

    % Define arguments for the smoothed_aggregation_solver constructor
    solver_constructor_args = {A_csr};

    % Convert B_nullspace to numpy array and add to constructor arguments if provided
    if ~isempty(B_nullspace_matlab)
        % Ensure B_nullspace_matlab is double for numpy conversion
        B_np_py = np.array(double(B_nullspace_matlab));
        % Add B to the solver constructor arguments as a keyword argument
        solver_constructor_args{end+1} = pyargs('B', B_np_py,'max_levels', int32(5));
        
        % --- Debugging: Display B_nullspace shape ---
        disp('Python B_nullspace shape:');
        disp(B_np_py.shape);
        % --- End Debugging ---
    end

%     solver_constructor_args{end+1} = pyargs('max_levels', 2);

    % --- Debugging: Display A_csr and b_np shapes before solving ---
    disp('Python A_csr shape:');
    disp(A_csr.shape);
    disp('Python b_np shape:');
    disp(b_np.shape);
    % --- End Debugging ---

    % Setup and solve using Smoothed Aggregation AMG solver
    % Pass A_csr and optionally B_nullspace to the solver constructor
    tic
    solver = py.pyamg.smoothed_aggregation_solver(solver_constructor_args{:}); 
    toc
    % Pass the hardcoded tolerance and maxiter, and the callback to the solve method
    tic
    x_py = solver.solve(b_np, solve_method_options);
    toc
    % Convert solution and residuals back to MATLAB arrays
    x = double(x_py);
    % Convert Python list of residuals to MATLAB array
    residualPCG = double(np.array(py.list(tracker.residuals)));
    save('/home/raul/Documents/GitHub/Article graphs/Por3D_AMG.mat','residualPCG')
%     residuals = double(py.list(tracker.residuals));

end
