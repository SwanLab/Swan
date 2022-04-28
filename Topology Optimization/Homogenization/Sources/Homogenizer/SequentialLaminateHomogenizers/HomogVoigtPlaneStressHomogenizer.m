classdef HomogVoigtPlaneStressHomogenizer < SequentialLaminateHomogenizer

    properties (Access = private)
    end

    methods (Access = public)

        function obj = HomogVoigtPlaneStressHomogenizer(c0,c1,dir,mi,theta)
            obj.init(c0,c1,dir,mi,theta);
            obj.computeAnisotropicTensors();
            obj.homogenize()
            obj.transformToVoigt()
            obj.makeHomogenizedTensorPlaneStress();
        end
    end

    methods (Access = protected)

        function transformToVoigt(obj)
            t = obj.homogTensor;
            t = Tensor2VoigtConverter.convert(t);
            obj.homogTensor = t;
        end

    end


end

