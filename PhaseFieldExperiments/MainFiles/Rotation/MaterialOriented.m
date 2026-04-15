classdef MaterialOriented < Material

    properties (Access = private)
        material
        mesh
        theta
    end

    methods (Access = public)

        function obj = MaterialOriented(cParams)
            obj.init(cParams);
        end

        function updateRotation(obj,thetaIn)
            obj.theta = thetaIn;
        end

        function C = obtainTensor(obj,phi)
            C0 = obj.material.obtainTensor(phi);
            s.operation = @(xV) obj.applyRotation(C0,xV);
            s.mesh = obj.mesh;
            s.ndimf = 16;
            C = DomainFunction(s);

            C.evaluate([0;0])
        end

        function dC = obtainTensorDerivative(obj,phi)
            dC0 = obj.material.obtainTensorDerivative(phi);
            s.operation = @(xV) obj.applyRotation(dC0,xV);
            s.mesh = obj.mesh;
            s.ndimf = 16;
            dC = DomainFunction(s);
        end

        function ddC = obtainTensorSecondDerivative(obj,phi)
            ddC0 = obj.material.obtainTensorSecondDerivative(phi);
            s.operation = @(xV) obj.applyRotation(ddC0,xV);
            s.mesh = obj.mesh;
            s.ndimf = 16;
            ddC = DomainFunction(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.material = cParams.material;
            obj.mesh = cParams.mesh;
            obj.theta = cParams.theta;
        end

    end

    methods (Access = private)

        function Crot = applyRotation(obj,C0,xV)
            C0 = C0.evaluate(xV);
            [Reps,Rsig] = obj.computeRotationMatrices(Expand(obj.theta,2).evaluate(xV));
            C0Voigt = convert2Voigt(C0,'Constitutive');
            CrotVoigt = pagemtimes(Rsig,pagemtimes(C0Voigt,Reps));
            Crot = convert2Tensor(CrotVoigt,'Constitutive');
        end

        function [Reps,Rsig] = computeRotationMatrices(~,theta)
            Reps = [(1+cos(2*theta))/2  (1-cos(2*theta))/2  -sin(2*theta)/2;
                    (1-cos(2*theta))/2  (1+cos(2*theta))/2   sin(2*theta)/2;
                        sin(2*theta)      -sin(2*theta)      cos(2*theta)  ];
            Rsig = pagetranspose(Reps);
        end

    end

end