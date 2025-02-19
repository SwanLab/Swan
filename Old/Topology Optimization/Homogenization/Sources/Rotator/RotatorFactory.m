classdef RotatorFactory < handle
    
    properties (Access = private)
        tensor
        angle
        dir
        rotator
    end
    
    methods (Access = public)
        
        function rotator = create(obj,tensor,angle,direction)
            obj.init(tensor,angle,direction)
            obj.createRotator()
            rotator = obj.getRotator();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,tensor,angle,normalDirection)
            obj.tensor = tensor;
            obj.angle = angle;
            obj.dir = normalDirection;
        end
        
        function createRotator(obj)
            a = obj.angle;
            d = obj.dir;
            rot = [];
            
            if obj.isFourthOrder()
                if obj.isVoigt()
                    if obj.is3D()
                        rot = FourthOrderVoigt3DRotator(a,d);
                    elseif obj.isPlaneStress()
                        rot = FourthOrderVoigtPlaneStressRotator(a,d);
                    end
                end
            elseif obj.isSecondOrder()
                if obj.isStress()
                    if obj.isVoigt()
                        if obj.is3D()
                            rot = StressVoigt3DRotator(a,d);
                        elseif obj.isPlaneStress()
                            rot = StressVoigtPlaneStressRotator(a,d);
                        end
                    elseif obj.isTensor()
                        if obj.is3D()
                            rot = StressRotator(a,d);
                        end
                    end
                end
                
            elseif obj.isVector()
                rot = VectorRotator(a,d);
            end
            
            if isempty(rot)
                error('Not admitted object to be Rotated')
            else
                obj.rotator = rot;
            end
        end
        
        function itIs = isVector(obj)
           itIs = strcmp(obj.tensor.getOrder(),'first');
        end
        
        function itIs = isStress(obj)
            itIs = strcmp(obj.tensor.getFieldName(),'stress');
        end
        
        function itIs = isFourthOrder(obj)
            itIs = strcmp(obj.tensor.getOrder(),'fourth');
        end
        
        function itIs = isSecondOrder(obj)
            itIs = strcmp(obj.tensor.getOrder(),'second');
        end
        
        function itIs = isTensor(obj)
            itIs = strcmp(obj.tensor.getRepresentation(), 'tensor');
        end
        
        function itIs = isVoigt(obj)
            itIs = strcmp(obj.tensor.getRepresentation(), 'voigt');
        end
        
        function itIs = is3D(obj)
            itIs = strcmp(obj.tensor.getElasticityCase(), '3D');
        end
        
        function itIs = isPlaneStress(obj)
            itIs = strcmp(obj.tensor.getElasticityCase(), 'planeStress');
        end
        
        function r = getRotator(obj)
            r = obj.rotator;
        end
        
    end
    
end

