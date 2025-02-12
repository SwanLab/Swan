classdef SymbolicVoigtRotationMatrixGenerator < handle
    
    properties (Access = public)
        VoigtRotationMatrix
    end
    
    
    properties (Access = private)
        RotationMatrix
        Stress
        RotatedStress
    end
    
    methods (Access = public)
        
        function obj = SymbolicVoigtRotationMatrixGenerator()
            obj.createRotationMatrix();
            obj.createStressTensor();
            obj.computeRotatedStressTensor();
            obj.obtainVoigtRotationMatrix()
        end
        
    end
    
    methods (Access = private)
        
        function obj = createRotationMatrix(obj)
            theta    = obj.createRotationAngle();
            vect     = obj.createNormalVector();            
            rotator = VectorRotator(theta,vect);
            obj.RotationMatrix = rotator.getRotationMatrix();
        end
        
        function u = createNormalVector(obj)
            u = sym('u',[3 1],'real');
        end
        
        function theta = createRotationAngle(obj)
            theta = sym('theta','real');
        end
        
        function createStressTensor(obj)
            obj.Stress = StressTensor();
            obj.Stress.tensor = sym('s',[3 3],'real');
            Tens = obj.Stress.tensor;
            Tens = obj.Stress.symmetrizeWithUpperDiagonal(Tens);
            obj.Stress.tensor = Tens;            
            obj.Stress.transformTensor2Voigt();
        end
        
        function computeRotatedStressTensor(obj)
            R = obj.RotationMatrix;
            S = obj.Stress.tensor;
            RotS = R'*S*R;
            obj.RotatedStress = StressTensor();
            obj.RotatedStress.tensor = simplify(RotS);
            obj.RotatedStress.transformTensor2Voigt();
        end
        
        function obtainVoigtRotationMatrix(obj)
            RotStre = obj.RotatedStress.tensorVoigt;
            Stre    = obj.Stress.tensorVoigt;
            DimVoigt = length(RotStre);
            RotMatrix = sym(zeros(DimVoigt,DimVoigt));
            for iStress = 1:DimVoigt
                RS = RotStre(iStress);
                for jStress = 1:DimVoigt
                    S = Stre(jStress);
                    [TensorValue,~] = coeffs(RS,S);
                    TensorValue = obj.takeApropiateComponent(TensorValue);
                    RotMatrix(iStress,jStress) = TensorValue;
                end
            end
            obj.VoigtRotationMatrix = simplify(RotMatrix);
        end
        
        function value = takeApropiateComponent(obj,value)
            if isempty(value)
                value = [];
            else
                Dim = length(value);
                if Dim == 1
                    value = 0;
                else
                    value = value(1);
                end
            end
        end
        
    end
    
end

