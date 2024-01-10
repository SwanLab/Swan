classdef VolumeConstraint < handle
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        mesh
        volumeTarget
        volume
    end
    
    methods (Access = public)        
        function  obj = VolumeConstraint(cParams)
            obj.init(cParams);
        end
        
        function computeFunctionAndGradient(obj,x)
            obj.volume.computeFunctionAndGradient(x);
            obj.computeFunction();
            obj.computeGradient();
        end  
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.volumeTarget = cParams.volumeTarget;
            obj.volume       = VolumeFunctional(cParams);
        end

        function computeFunction(obj)
            v         = obj.volume.value;
            vTar      = obj.volumeTarget;
            obj.value = v/vTar-1;
        end

        function computeGradient(obj)
            vTar         = obj.volumeTarget;
            vGrad        = obj.volume.gradient;
            s.mesh       = obj.mesh;
            s.feFunType  = class(vGrad);
            s.ndimf      = 1;
            g            = FeFunction.createEmpty(s);
            g.fValues    = vGrad.fValues/vTar;
            obj.gradient = g;
        end
    end

end

