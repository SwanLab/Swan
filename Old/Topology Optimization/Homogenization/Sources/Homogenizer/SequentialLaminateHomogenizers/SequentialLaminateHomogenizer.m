classdef SequentialLaminateHomogenizer < handle
    
     
    properties (Access = protected)
        directions
        laminateParams
        rank
        theta
        fiberTensor
        matrixTensor
        anisotropicTensors
        
        aniTensors
        planeStressHomogenizedTensor
        
        homogTensor
        fibTen
        matTen
    end
    
    methods (Access = public)

        function PS = getPlaneStressHomogenizedTensor(obj)
            PS = obj.planeStressHomogenizedTensor;
        end 
        
    end
    
    methods (Access = protected)
        function init(obj,C0,C1,dir,mi,theta)
            obj.fibTen         = C1;
            obj.matTen         = C0;
            obj.laminateParams = mi;
            obj.rank           = length(obj.laminateParams);
            obj.theta          = theta;
            obj.storeDirections(dir)
            obj.checkConsistencyLaminateParamsAndDirections()
        end
        
        function storeDirections(obj,dir)
            for irank = 1:obj.rank
                d = dir{irank};
                d.normalize();
                obj.directions(irank,:) = d.getValue();
            end
        end
        
        function checkConsistencyLaminateParamsAndDirections(obj)
            assert(size(obj.directions,1) == length(obj.laminateParams))
        end
        
        function computeAnisotropicTensors(obj)
            for irank = 1:obj.rank
                c1 = obj.fibTen;
                dir = obj.directions(irank,:);
                obj.anisotropicTensors{irank} = AnisotropicContributionTensor(c1,dir);
            end
        end
        

        function homogenize(obj)
            c1 = obj.fibTen;
            c0 = obj.matTen;
            ani = obj.anisotropicTensors;
            mi = obj.laminateParams;
            angle = obj.theta;
            hom = Homogenizer(c0,c1,ani,mi,angle);
            obj.homogTensor = hom.getHomogenizedTensor();
        end
        
         function makeHomogenizedTensorPlaneStress(obj)
            ht = obj.homogTensor;
            ps = PlaneStressTransformer.transform(ht);
            obj.planeStressHomogenizedTensor = ps;
        end
        
         function transformTensorsToVoigt(obj)
            obj.transforMatrixAndFiberTensorInVoigt();
            obj.transformAnisotropicTensorsInVoigtNotation();
         end
         
    end
    
    methods (Access = private)
        
        function transforMatrixAndFiberTensorInVoigt(obj)
            obj.matTen = Tensor2VoigtConverter.convert(obj.matTen);
            obj.fibTen = Tensor2VoigtConverter.convert(obj.fibTen);
        end
        
        function transformAnisotropicTensorsInVoigtNotation(obj)
            for irank = 1:obj.rank
                ani  = obj.anisotropicTensors{irank};
                vAni = Tensor2VoigtConverter.convert(ani);
                obj.anisotropicTensors{irank} = vAni;
            end

        end
    end
    
end

