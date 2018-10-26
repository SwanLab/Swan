classdef TensorInverter < Inverter
    
    properties (Access = protected)
        invertedTensor
    end
    
    properties (Access = private)
       dim 
       tensor            
       matrix
       invMatrix
    end
    
    methods (Access = public)
        
        function obj = TensorInverter(tensor)
            obj.init(tensor)
            obj.transformTensor2Matrix()
            obj.makeInverseOfMatrix()
            obj.transformMatrix2Tensor()
        end        
    end
    
    methods (Access = private)
        
        function init(obj,tensor)
           obj.tensor = tensor;             
           obj.dim = size(obj.tensor,1);
        end
        
        function transformTensor2Matrix(obj)
            d = obj.dim;
            a = obj.tensor;
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
           d = obj.dim; 
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
           obj.invertedTensor = a;
        end
        
    end
end

