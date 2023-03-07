classdef StifnessMatrixComputer < handle
    properties (Access = public)
        globalStifnessMatrix
    end
    properties (Access = private)
%        elementType
        t
        poissonCoefficient
        elasticModuleMinimun
        elasticModuleNeutral
        projectedField
        penalization
        elementNumberX
        elementNumberY
        conectivityMatrixMat

        elementalStiffnessMatrix
        sensitizedElasticModule
    end
    methods (Access = public)
        function obj = StifnessMatrixComputer(cParams)
            obj.inputData(cParams)
        end

        function compute(obj)
%            obj.computeElementalStiffnessMatrices()
            obj.penalizeElasticModule()
            obj.computeGlobalStifnessMatrix()
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.elementalStiffnessMatrix = cParams.structure.elementalStiffnessMatrix;
            obj.t=cParams.structure.t;
            obj.poissonCoefficient=cParams.structure.poissonCoefficient;
            obj.elasticModuleMinimun=cParams.structure.elasticModuleMinimun;
            obj.elasticModuleNeutral=cParams.structure.elasticModuleNeutral;  
            obj.penalization=cParams.structure.penalization;
            
            obj.elementNumberX = cParams.mesh.elementNumberX;
            obj.elementNumberY = cParams.mesh.elementNumberY;
            obj.conectivityMatrixMat = cParams.mesh.conectivityMatrixMat;

            obj.projectedField=cParams.projectedField;
        end
%         function computeElementalStiffnessMatrices(obj)
%             s.elementType = obj.elementType;
%             s.t = obj.t;
%             s.poissonCoefficient = obj.poissonCoefficient;
%             B = ElementalStiffnessMatricesComputer(s);
%             B.compute();
%             obj.elementalStiffnessMatrix = B.elementalStiffnessMatrix;
%         end
        function penalizeElasticModule(obj)
            s.elasticModuleMinimun = obj.elasticModuleMinimun;
            s.elasticModuleNeutral = obj.elasticModuleNeutral;
            s.projectedField = obj.projectedField;
            s.penalization = obj.penalization;
            s.nonPenalizedVariable = obj.elementalStiffnessMatrix;
            B = Penalizer(s);
            B.penalize();
            obj.elementalStiffnessMatrix = B.penalizedVariable;
        end
        function computeGlobalStifnessMatrix(obj)
            iK      = reshape(kron(obj.conectivityMatrixMat,ones(8,1))',64*obj.elementNumberX*obj.elementNumberY,1);
            jK      = reshape(kron(obj.conectivityMatrixMat,ones(1,8))',64*obj.elementNumberX*obj.elementNumberY,1);
            sKiI  = reshape(obj.elementalStiffnessMatrix,64*obj.elementNumberX*obj.elementNumberY,1);
            obj.globalStifnessMatrix   = sparse(iK,jK,sKiI); 
            obj.globalStifnessMatrix = (obj.globalStifnessMatrix+obj.globalStifnessMatrix')/2;
        end 
    end
end