classdef RHSIntegratorFactory < handle

    methods (Access = public, Static)

        function obj = create(cParams)
            switch cParams.type
                case 'ShapeFunction'
                    obj = RHSIntegratorShapeFunction(cParams);
                case 'ShapeFunctionRT'
                    obj = RHSintegrator_ShapeFunctionRT(cParams);
                case 'ShapeFunctionN'
                    obj = RHSintegrator_ShapeFunctionN(cParams);
                case 'ShapeDerivative'
                    obj = RHSIntegratorShapeDerivative(cParams);
                case 'ShapeDerivativeTensor'
                    obj = RHSIntegratorShapeDerivativeTensor(cParams);
                case 'ShapeSymmetricDerivative'
                    obj = RHSIntegratorShapeSymmDerivative(cParams);
                case 'CutMesh'
                    obj = RHSintegrator_CutMesh(cParams);
                case {'Composite','COMPOSITE'}
                    obj = RHSintegrator_Composite(cParams);
                case 'Elastic'
                    switch cParams.scale
                        case 'MACRO'
                            obj = RHSIntegratorElasticMacro(cParams);
                        case 'MICRO'
                            obj = RHSIntegratorElasticMicro(cParams);
                    end
                case 'ElasticMicro'
                    obj = RHSIntegratorElasticMicro(cParams);
                case 'Stokes'
                    obj = RHSIntegratorStokes(cParams);
                case 'Unfitted'
                    obj = RHSIntegratorUnfitted(cParams);
            end
        end

    end

end