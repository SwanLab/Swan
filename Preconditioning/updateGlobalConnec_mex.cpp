#include "mex.hpp"
#include "mexAdapter.hpp"
#include <vector>
#include <algorithm>

using namespace matlab::mex;
using namespace matlab::data;

class MexFunction : public matlab::mex::Function {
    public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {

        // Check input/output counts
        if (inputs.size() != 2) {
            getEngine()->feval(u"error", 0,
                               std::vector<Array>({ factory.createScalar(u"updateGlobalConnec_mex: Two inputs required: connecGlob and interfaceConnec.") }));
        }
        if (outputs.size() > 1) {
            getEngine()->feval(u"error", 0,
                               std::vector<Array>({ factory.createScalar(u"updateGlobalConnec_mex: Only one output allowed.") }));
        }

        // Get input arrays (double from MATLAB)
        TypedArray<double> connecGlob = inputs[0];
        TypedArray<double> interfaceConnec = inputs[1];

        // Get dimensions
        auto dims = connecGlob.getDimensions();
        size_t nref = interfaceConnec.getDimensions()[0];
        size_t ncopy = interfaceConnec.getDimensions()[1];

        // Convert inputs to int32 for faster processing
        std::vector<int32_t> updtConnecGlob(connecGlob.begin(), connecGlob.end());
        std::vector<int32_t> rC(interfaceConnec.begin(), interfaceConnec.end());

        // Main algorithm
        for (size_t iref = 0; iref < nref; ++iref) {
            // Collect non-zero entries from the row
            std::vector<int32_t> aux;
            for (size_t icopy = 0; icopy < ncopy; ++icopy) {
                int32_t val = rC[iref + icopy * nref]; // MATLAB column-major indexing
                if (val > 0) aux.push_back(val);
            }

            for (size_t icopy = 1; icopy < aux.size(); ++icopy) {
                int32_t target = aux[icopy];
                int32_t first = aux[0];

                // Update updtConnecGlob
                for (auto &v : updtConnecGlob) {
                    if (v == target) v = first;
                    else if (v > target) v -= 1;
                }

                // Update aux
                for (auto &v : aux) {
                    if (v > target) v -= 1;
                }

                // Update rC
                for (auto &v : rC) {
                    if (v > target) v -= 1;
                }
            }
        }


        std::vector<double> updtConnecGlobDouble(updtConnecGlob.begin(), updtConnecGlob.end());

        // Create MATLAB array from double vector
        TypedArray<double> out = factory.createArray(dims, updtConnecGlobDouble.begin(), updtConnecGlobDouble.end());
        outputs[0] = out;
    }

    private:
    ArrayFactory factory;
};
