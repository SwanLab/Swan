// Call Function generated automatically on 27-May-2025 01:35:29
// Created on DESKTOP-BM120QM
#include "cstdint"
#include "mex.h"
#include "math.h"
#include "algorithm"
#include "array"
#include "mb_fm_mex_source_zermelo.h"
#include "mb_fm_mex_source_zermelo_types.h"

// Constants - Input Dimensions
#define DIM_M_STATES 2
#define DIM_M_CONTROLS 1

// Constants - Number of Independent Variables (non discrete control case)
#define N_IDP 3

// Constants - Output Sizes
#define NUM_OUT_STATESDOT 2

static int const nthreads = 1;

// Function Header
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) try {

    // Declare input arrays
    mxArray const *array_states = prhs[0];
    mxArray const *array_controls = prhs[1];

    // Declare number of time steps
    int nEval;

    // Declare Outputs
    double *statesdot;
    double *j_statesdot;
    mwSize j_dim_statesdot[3];

    // Write Call Information
    if (nrhs == 0 && nlhs == 0) {
        // System Information
        mexPrintf("<strong>System Information</strong>\n");
        mexPrintf("Mex file created by Derivative Model Builder\n");
        mexPrintf("- Date:                27-May-2025 01:35:29\n");
        mexPrintf("- Computer:            DESKTOP-BM120QM\n");
        mexPrintf("- MATLAB-Version:      9.14.0.2891782 (R2023a) Update 8\n");
        mexPrintf("- DerivativeOrder:     1\n");
        mexPrintf("- Jacobian Calculation:true\n");
        mexPrintf("- Hessian-Calculation: false\n");
        mexPrintf("\n");

        // Input Information
        mexPrintf("<strong>Input Information</strong>\n");
        mexPrintf("<strong>Name</strong>     <strong>Size</strong>  <strong>DataType</strong> <strong>Derivative</strong> <strong>MultipleTimeEval</strong> <strong>VariableSize</strong>\n");
        mexPrintf("states   [2 1] double       true             true        false\n");
        mexPrintf("| x\n");
        mexPrintf("| y\n");
        mexPrintf("\n");
        mexPrintf("controls [1 1] double       true             true        false\n");
        mexPrintf("| u\n");
        mexPrintf("\n");
        mexPrintf("\n");

        // Output Information
        mexPrintf("<strong>Output Information</strong>\n");
        mexPrintf("<strong>Name</strong>      <strong>Size</strong> \n");
        mexPrintf("statesdot [2 1]\n");

        return;
    }
    else if (nrhs == 0 && nlhs == 1) {
        const char *names[] = {"input", "output", "info", "name", "type", "WrapperClass"};
        const char *i_names[] = {"m", "n","name","argnames","type", "groupindex", "DataType"};
        const char *o_names[] = {"m", "n","name","argnames","type","jac_sparsity","hess_sparsity"};
        mxArray *struct_inputs;
        mxArray *struct_outputs;
        const char *info_names[] = {"Date", "Computer", "MATLAB", "DerivativeOrder", "Jacobian", "Hessian", "UseSparsityEstimator"};
        mxArray *struct_info;
        mxArray *mx;
        double *sparsity_j;
        double *sparsity_h;

        plhs[0] = mxCreateStructMatrix(1,1,6,names);
        mxSetField(plhs[0], 0, names[3], mxCreateString("fm_mex_source_zermelo"));
        mxSetField(plhs[0], 0, names[4], mxCreateString("SIMULATION_MODEL"));
        mxSetField(plhs[0], 0, names[5], mxCreateString("falcon.ModelWrapper"));

        struct_inputs = mxCreateStructMatrix(2,1,7,i_names);

        mxSetField(struct_inputs, 0, i_names[0], mxCreateDoubleScalar(2));
        mxSetField(struct_inputs, 0, i_names[1], mxCreateDoubleScalar(1));
        mxSetField(struct_inputs, 0, i_names[2], mxCreateString("states"));
        mxSetField(struct_inputs, 0, i_names[4], mxCreateString("STATE"));
        mxSetField(struct_inputs, 0, i_names[5], mxCreateDoubleScalar(0));
        mxSetField(struct_inputs, 0, i_names[6], mxCreateString("double"));
        mx = mxCreateCellMatrix(2, 1);
        mxSetCell(mx, 0,  mxCreateString("x"));
        mxSetCell(mx, 1,  mxCreateString("y"));
        mxSetField(struct_inputs, 0, i_names[3], mx);

        mxSetField(struct_inputs, 1, i_names[0], mxCreateDoubleScalar(1));
        mxSetField(struct_inputs, 1, i_names[1], mxCreateDoubleScalar(1));
        mxSetField(struct_inputs, 1, i_names[2], mxCreateString("controls"));
        mxSetField(struct_inputs, 1, i_names[4], mxCreateString("CONTROL"));
        mxSetField(struct_inputs, 1, i_names[5], mxCreateDoubleScalar(0));
        mxSetField(struct_inputs, 1, i_names[6], mxCreateString("double"));
        mx = mxCreateCellMatrix(1, 1);
        mxSetCell(mx, 0,  mxCreateString("u"));
        mxSetField(struct_inputs, 1, i_names[3], mx);
        mxSetField(plhs[0], 0, names[0], struct_inputs);

        struct_outputs = mxCreateStructMatrix(1,1,7,o_names);

        mxSetField(struct_outputs, 0, o_names[0], mxCreateDoubleScalar(2));
        mxSetField(struct_outputs, 0, o_names[1], mxCreateDoubleScalar(1));
        mxSetField(struct_outputs, 0, o_names[2], mxCreateString("statesdot"));
        mxSetField(struct_outputs, 0, o_names[4], mxCreateString("STATEDOT"));
        mx = mxCreateDoubleMatrix(2, 3, mxREAL);
        sparsity_j = reinterpret_cast<double *>(mxGetPr(mx));
        {
            std::array<bool, 2ul * 3ul> const sparsity_j_template = {
                0, 0,  /* column 1 */
                1, 0,  /* column 2 */
                1, 1,  /* column 3 */
            };  /* end of sparsity_j_template */
            std::copy(sparsity_j_template.cbegin(), sparsity_j_template.cend(), sparsity_j);
        }
        mxSetField(struct_outputs, 0, o_names[5], mx);
        mx = mxCreateDoubleMatrix(0, 0, mxREAL);
        sparsity_h = reinterpret_cast<double *>(mxGetPr(mx));
        {
            std::array<bool, 0ul * 0ul> const sparsity_h_template = {
            };  /* end of sparsity_h_template */
            std::copy(sparsity_h_template.cbegin(), sparsity_h_template.cend(), sparsity_h);
        }
        mxSetField(struct_outputs, 0, o_names[6], mx);
        mx = mxCreateCellMatrix(0,0);
        mxSetField(struct_outputs, 0, o_names[3], mx);
        mxSetField(plhs[0], 0, names[1], struct_outputs);
        struct_info = mxCreateStructMatrix(1,1,7,info_names);
        mxSetField(struct_info, 0, info_names[0], mxCreateString("27-May-2025 01:35:29"));
        mxSetField(struct_info, 0, info_names[1], mxCreateString("DESKTOP-BM120QM"));
        mxSetField(struct_info, 0, info_names[2], mxCreateString("9.14.0.2891782 (R2023a) Update 8"));
        mxSetField(struct_info, 0, info_names[3], mxCreateDoubleScalar(1));
        mxSetField(struct_info, 0, info_names[4], mxCreateLogicalScalar(true));
        mxSetField(struct_info, 0, info_names[5], mxCreateLogicalScalar(false));
        mxSetField(struct_info, 0, info_names[6], mxCreateLogicalScalar(false));
        mxSetField(plhs[0], 0, names[2], struct_info);
        return;
    }

    // Check Number of Input Arguments
    if (nrhs < 2 || nrhs > 2) {
        mexErrMsgIdAndTxt("MATLAB:callModel:Error","Wrong number of input arguments for fm_mex_source_zermelo. The required number of input arguments is 2. Call the mex file with no input arguments to get exact input & output type information.");
    }

    // Extract Inputs

    if (!mxIsClass(array_states, "double") || mxIsComplex(array_states)) {
        mexErrMsgIdAndTxt("MATLAB:callModel:Error", "Wrong argument type: states must be real double, found %s.", mxGetClassName(array_states));
    }
    double const *states = reinterpret_cast<double const *>(mxGetPr(array_states));

    if (!mxIsClass(array_controls, "double") || mxIsComplex(array_controls)) {
        mexErrMsgIdAndTxt("MATLAB:callModel:Error", "Wrong argument type: controls must be real double, found %s.", mxGetClassName(array_controls));
    }
    double const *controls = reinterpret_cast<double const *>(mxGetPr(array_controls));
    nEval = mxGetN(array_states);

    // Check Input Dimensions
    // Input_1: states
    if (mxGetM(array_states) != DIM_M_STATES) {
        mexErrMsgIdAndTxt("MATLAB:callModel:Error", "Number of rows for input 1 with name \"states\" should be 2 but is %i.", mxGetM(prhs[0]));
    }
    if (mxGetN(array_states) != nEval) {
        mexErrMsgIdAndTxt("MATLAB:callModel:Error", "Number of columns for input 1 with name \"states\" should be equal to the number of time steps (%i) but is %i.", nEval, mxGetN(prhs[0]));
    }
    // Input_2: controls
    if (mxGetM(array_controls) != DIM_M_CONTROLS) {
        mexErrMsgIdAndTxt("MATLAB:callModel:Error", "Number of rows for input 2 with name \"controls\" should be 1 but is %i.", mxGetM(prhs[1]));
    }
    if (mxGetN(array_controls) != nEval) {
        mexErrMsgIdAndTxt("MATLAB:callModel:Error", "Number of columns for input 2 with name \"controls\" should be equal to the number of time steps (%i) but is %i.", nEval, mxGetN(prhs[1]));
    }

    // Prepare jacobian and hessian dimensions
    // Output_1: statesdot
    plhs[0] = mxCreateDoubleMatrix(2, 1*nEval, mxREAL);
    statesdot = mxGetPr(plhs[0]);
    // Output_jacobian_1: statesdot
    j_dim_statesdot[0] = NUM_OUT_STATESDOT;
    j_dim_statesdot[1] = N_IDP;
    j_dim_statesdot[2] = nEval;
    plhs[1] = mxCreateNumericArray(3,j_dim_statesdot,mxDOUBLE_CLASS,mxREAL);
    j_statesdot = mxGetPr(plhs[1]);

    // Call Model in for-loop
    for (int iEval = 0; iEval < nEval; iEval++) {
        mb_fm_mex_source_zermelo(
        &states[DIM_M_STATES * iEval],
        controls[DIM_M_CONTROLS * iEval],
        &statesdot[NUM_OUT_STATESDOT * iEval],
        &j_statesdot[NUM_OUT_STATESDOT * N_IDP * iEval]
        );
    }
    // Remove Temporary Variables - Variable Size Data

}
catch (std::exception& e)
{
    mexErrMsgIdAndTxt("MATLAB:callModel:Error", "%s", e.what());
}
