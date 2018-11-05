classdef TensorInverter < Inverter
    
    properties (Access = private)
       matrix
       invMatrix
    end
    
    methods (Access = public)
        
        function obj = TensorInverter(tensor)
            obj.compute(tensor);
        end        
    end
    
    methods (Access = protected)
        
        function computeInverse(obj)
            obj.transformTensor2Matrix()
            obj.makeInverseOfMatrix()
            obj.transformMatrix2Tensor()
        end
    end
    
    
    methods (Access = private)
        
        function transformTensor2Matrix(obj)
            d = obj.tensor.getDimension();
            a = obj.tensor.getValue();
            am = zeros(d*d,d*d);
            for i = 1:d
                for j = 1:d
                    for k = 1:d
                        for l = 1:d
                            im = d*(i - 1) + j;
                            jm = d*(k - 1) + l;
                            am(im,jm)  = a(i,j,k,l);
                        end
                    end
                end
            end
            obj.matrix = am;
        end
        
        function makeInverseOfMatrix(obj)
            obj.invMatrix = inv(obj.matrix);            
        end
           
        function transformMatrix2Tensor(obj)
           d = obj.tensor.getDimension();
           am = obj.invMatrix;
           a = zeros(d,d,d,d);           
           for i = 1:d
               for j = 1:d
                   for k = 1:d
                       for l = 1:d
                           im = d*(i - 1) + j;    
                           jm = d*(k - 1) + l;
                           a(i,j,k,l) = am(im,jm);
                       end
                   end
               end
           end           
           obj.invertedTensor.setValue(a);
        end
        
    end
end

