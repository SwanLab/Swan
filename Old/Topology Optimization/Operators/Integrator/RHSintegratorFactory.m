classdef RHSintegratorFactory < handle

    methods (Access = public, Static)

        function obj = create(cParams)
            switch cParams.type
                case 'ShapeFunction'
                    % Computes the RHS using the NODAL FORCES and SHAPE
                    % FUNCTIONS
                    obj = RHSintegrator_ShapeFunction(cParams);
                case 'ShapeFunctionRT'
                    % Computes the RHS using the NODAL FORCES and SHAPE
                    % FUNCTIONS
                    obj = RHSintegrator_ShapeFunctionRT(cParams);
                case 'ShapeFunctionN'
                    % Computes the RHS using the NODAL FORCES and SHAPE
                    % FUNCTIONS
                    obj = RHSintegrator_ShapeFunctionN(cParams);
                case 'ShapeDerivative'
                    % Computes the RHS using the NODAL FORCES and SHAPE
                    % FUNCTIONS' DERIVATIVE
                    obj = RHSintegrator_ShapeDerivative(cParams);
                case 'ShapeSymmetricDerivative'
                    % Computes the RHS using the NODAL FORCES and SHAPE
                    % FUNCIONS' SYMMETRIC GRADIENT
                    obj = RHSintegrator_ShapeSymmDerivative(cParams);
                case 'CutMesh'
                    % Computes the RHS using the NODAL FORCES and SHAPE
                    % FUNCTIONS for CUT meshes
                    obj = RHSintegrator_CutMesh(cParams);
                case {'Composite','COMPOSITE'}
                    % Computes the RHS as a COMPOSITE of other RHS
                    % integrators
                    obj = RHSintegrator_Composite(cParams);
                    
                case 'Elastic'
                    % Computes the RHS for ELASTIC problems
                    switch cParams.scale
                        case 'MACRO'
                            obj = RHSintegrator_ElasticMacro(cParams);
                        case 'MICRO'
                            obj = RHSintegrator_ElasticMicro(cParams);
                    end
                case 'ElasticMicro'
                    obj = RHSintegrator_ElasticMicro(cParams);
                case 'Stokes'
                    obj = RHSintegrator_Stokes(cParams);
                case 'Unfitted'
                    obj = RHSIntegratorUnfitted(cParams);
            end
        end

    end

end