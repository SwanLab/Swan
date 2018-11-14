classdef testExplicitAndSymbolicVoigtRotation < test
    
    properties (Access = private)
       Rotation
       SymRotation
    end
    
    methods (Access = public)
        
        function obj = testExplicitAndSymbolicVoigtRotation()
            obj.createSymVoigtRotationMatrix();
            obj.createVoigtRotationMatrix()    ;                       
        end
        
    end
    
    methods (Access = private)
        
        function createSymVoigtRotationMatrix(obj)
           SymMatrixGenerator = SymbolicVoigtRotationMatrixGenerator();
           obj.SymRotation = SymMatrixGenerator.VoigtRotationMatrix();            
        end
        
        function createVoigtRotationMatrix(obj)
           Angle = sym('theta','real');
           Direction = sym('u',[3 1],'real');
           MatrixGenerator = VoigtRotationMatrixGenerator(Angle,Direction);
           obj.Rotation = MatrixGenerator.VoigtMatrix();
        end
    end
    
    
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)   
            hasPassed = double(norm(obj.Rotation - obj.SymRotation) )< 1e-14;
        end
        
    end
end

