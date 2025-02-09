classdef Inverter < handle
    
    properties (Access = protected)
        invertedTensor
        tensor
    end
    
    methods (Access = public,Static)
        
        function invertedTensor = invert(Tensor)
            factory  = InverterFactory();
            inverter = factory.create(Tensor);
            invertedTensor = inverter.getInvertedTensor();
        end
    end
    
    methods (Access = protected)
        
        function compute(obj,tensor)
            obj.init(tensor)
            obj.createInvertedTensor()
            obj.computeInverse()
        end
        
    end
    
    methods (Access = private)

        function createInvertedTensor(obj)
            if isa(obj.tensor,'FourthOrder3DTensor')
                obj.invertedTensor = obj.tensor.clone();
            else
                
            if strcmp(obj.tensor.getFieldName,'stiffness')
                if  strcmp(obj.tensor.getElasticityCase,'3D')
                    if  strcmp(obj.tensor.getRepresentation(),'voigt')
                        obj.invertedTensor = Compliance3DVoigtTensor();
                    else
                        obj.invertedTensor = Compliance3DTensor();
                    end
                elseif strcmp(obj.tensor.getElasticityCase,'planeStress')
                    if strcmp(obj.tensor.getRepresentation(),'voigt')
                        obj.invertedTensor = CompliancePlaneStressVoigtTensor();
                    else
                        obj.invertedTensor = CompliancePlaneStressTensor();
                    end
                end
            elseif strcmp(obj.tensor.getFieldName,'compliance')
                if  strcmp(obj.tensor.getElasticityCase,'3D')
                    if strcmp(obj.tensor.getRepresentation(),'voigt')
                       obj.invertedTensor = Stiffness3DVoigtTensor();
                    else
                       obj.invertedTensor = Stiffness3DTensor();
                    end
                elseif strcmp(obj.tensor.getElasticityCase,'planeStress')
                    if strcmp(obj.tensor.getRepresentation(),'voigt')
                       obj.invertedTensor = StiffnessPlaneStressVoigtTensor();
                    else
                       obj.invertedTensor = StiffnessPlaneStressTensor();
                    end
                end
            else
                obj.invertedTensor = obj.tensor.clone();
            end
            end
        end
        
        function T = getInvertedTensor(obj)
            T = obj.invertedTensor;
        end
        
        function init(obj,t)
            obj.tensor = t;
        end
        
    end
    
    methods (Abstract, Access = protected)
       computeInverse(obj)
    end

end

