classdef ElementalStiffnessMatricesComputer < handle

    properties (Access = public)
        elementalStiffnessMatrix
    end
    properties (Access = private)
        poissonCoefficient
        t
        elementType
    end
    methods (Access = public)
        function obj = ElementalStiffnessMatricesComputer(cParams)
            obj.inputData(cParams)
        end
        function compute(obj)
            obj.chose()
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.elementType = cParams.mesh.elementType;
            obj.t = cParams.t;
            obj.poissonCoefficient = cParams.poissonCoefficient;
        end
        function chose(obj)
            switch obj.elementType
                case 'QUAD'
                    obj.computeSquareElement();
                case 'TRIANGLE'
                    obj.computeTriangularElement();
                otherwise
                    disp('Element mesh not implemented')
            end
        end
        function computeSquareElement(obj)
            s.poissonCoefficient = obj.poissonCoefficient;
            s.t = obj.t;
            B = SquareElement(s);
            B.computeStifnessMatrix();
            obj.elementalStiffnessMatrix = B.stifnessMatrix;     
        end
        function computeTriangularElement(obj)
            s.poissonCoefficient = obj.poissonCoefficient;
            s.t = obj.t;
            B = TriangularElement(s);
            B.computeStifnessMatrix();
            obj.elementalStiffnessMatrix = B.stifnessMatrix;     
        end
       
    end
end