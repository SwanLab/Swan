classdef ProductComputerOutputCreator < handle
    
    properties (Access = private)
        fourthOrder
        secondOrder
        secondOrderOut
    end
    
    methods (Access = public)
        
        function obj = ProductComputerOutputCreator(C,e)
            obj.init(C,e)
            obj.createSecondOrderOut()
        end
        
        function s = getSecondOrderOut(obj)
            s = obj.secondOrderOut;
        end
    end
    
    
    methods (Access = private)
        
        function init(obj,C,e)
            obj.fourthOrder = C;
            obj.secondOrder = e;            
        end
        
        function createSecondOrderOut(obj)
            
            if obj.isTensor()
                if obj.is3D()
                    if obj.isStress()
                        s = Stress3DTensor;
                    elseif obj.isStrain()
                        s = Strain3DTensor;
                    end
                elseif obj.isPlaneStress()
                    if obj.isStress()
                        s = StressPlaneStressTensor;
                    elseif obj.isStrain()
                        s = StrainPlaneStressTensor;
                    end
                end
            elseif obj.isVoigt()
                if obj.is3D()
                    if obj.isStress()
                        s = Stress3DVoigtTensor();
                    elseif obj.isStrain()
                        s = Strain3DTensor;
                    end
                elseif obj.isPlaneStress()
                    if obj.isStress()
                        s = StressPlaneStressVoigtTensor;
                    elseif obj.isStrain()
                        s = StrainPlaneStressVoigtTensor;
                    end
                end
            end
            obj.secondOrderOut = s;

        end
        
        function itIs = isVoigt(obj)
            isFourthVoigt = strcmp(obj.fourthOrder.getRepresentation(),'voigt');
            isSecondVoigt = strcmp(obj.secondOrder.getRepresentation(),'voigt');
            itIs = isFourthVoigt && isSecondVoigt;
        end
        
        function itIs = isTensor(obj)
            isFourthVoigt = strcmp(obj.fourthOrder.getRepresentation(),'tensor');
            isSecondVoigt = strcmp(obj.secondOrder.getRepresentation(),'tensor');
            itIs = isFourthVoigt && isSecondVoigt;
        end
        
        function itIs = isStress(obj)
            isConstiutive = strcmp(obj.fourthOrder.getOrder(),'fourth');
            isStrian = strcmp(obj.secondOrder.getFieldName(),'strain');
            itIs = isConstiutive && isStrian;                        
        end        
        
        function itIs = isStrain(obj)
            isConstiutive = strcmp(obj.fourthOrder.getOrder(),'fourth');
            isStress = strcmp(obj.secondOrder.getFieldName(),'stress');
            itIs = isConstiutive && isStress;
        end
        
        function itIs = is3D(obj)
            isFourth3D = strcmp(obj.fourthOrder.getElasticityCase(),'3D');
            isSecond3D = strcmp(obj.secondOrder.getElasticityCase(),'3D');
            itIs = isFourth3D && isSecond3D;
        end
        
        function itIs = isPlaneStress(obj)
            isFourthPS = strcmp(obj.fourthOrder.getElasticityCase(),'planeStress');
            isSecondPS = strcmp(obj.secondOrder.getElasticityCase(),'planeStress');
            itIs = isFourthPS && isSecondPS;
        end

        
    end
    
end

