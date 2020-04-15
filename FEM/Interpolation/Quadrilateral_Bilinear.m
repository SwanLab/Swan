classdef Quadrilateral_Bilinear < Interpolation
    
    methods (Access = public)
        
        function obj = Quadrilateral_Bilinear(cParams)
            obj.init(cParams);
            obj.computeParameters();
            obj.computeCases();
        end
        
        function computeShapeDeriv(obj,posgp)
            obj.computeShapes(posgp);
            obj.computeShapeDerivatives(posgp)
        end
        
    end
    
    methods (Access = private)
        
        function computeParameters(obj)
            obj.type = 'QUADRILATERAL';
            obj.ndime = 2;
            obj.nnode = 4;
            obj.pos_nodes = [-1 -1; 1 -1; 1 1; -1 1];
            obj.isoDv = 4;                        
        end
        
        function computeShapes(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);            
            s = posgp(1,:,:);
            t = posgp(2,:,:);
            I = ones(size(t));          
            obj.shape = zeros(obj.nnode,ngaus,nelem);
            obj.shape(1,:,:) = 0.25*(I-t-s+s.*t);
            obj.shape(2,:,:) = 0.25*(I-t+s-s.*t);
            obj.shape(3,:,:) = 0.25*(I+t+s+s.*t);
            obj.shape(4,:,:) = 0.25*(I+t-s-s.*t); 
        end
        
        function computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:);
            t = posgp(2,:,:);
            I = ones(size(t));
            obj.deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            obj.deriv(1,1,:,:) = 0.25*(-I+t);
            obj.deriv(1,2,:,:) = 0.25*(+I-t);
            obj.deriv(1,3,:,:) = 0.25*(+I+t);
            obj.deriv(1,4,:,:) = 0.25*(-I-t);
            obj.deriv(2,1,:,:) = 0.25*(-I+s);
            obj.deriv(2,2,:,:) = 0.25*(-I-s);
            obj.deriv(2,3,:,:) = 0.25*(+I+s);
            obj.deriv(2,4,:,:) = 0.25*(+I-s);            
        end
        
        function computeCases(obj)
            obj.iteration = [1 2 3 4;
                             2 3 4 1];
            obj.cases(:,:,1)=[1 5 6
                5 2 3
                5 3 6
                6 3 4;
                zeros(2,3)];
            obj.cases(:,:,2)=[1 5 4
                5 2 6
                5 6 4
                6 3 4;
                zeros(2,3)];
            obj.cases(:,:,3)=[1 2 5
                1 5 6
                5 3 6
                1 6 4;
                zeros(2,3)];
            obj.cases(:,:,4)=[1 2 6
                6 5 4
                2 3 5
                6 2 5;
                zeros(2,3)];
            obj.cases(:,:,5)=[1 5 4
                5 6 4
                5 2 6
                2 3 6;
                zeros(2,3)];
            obj.cases(:,:,6)=[6 3 4
                5 3 6
                1 5 6
                1 2 5;
                zeros(2,3)];
            obj.cases(:,:,7)=[1 5 8
                5 2 8
                8 2 4
                2 6 4
                4 6 7
                6 3 7];
            obj.main_loop=[4 3];
            obj.extra_cases=[7];
            %             obj.main_loop=[3 4]
            %             obj.extra_loop=[7]
            %             obj.cases{1}=[1 5 6
            %                 5 2 3
            %                 5 3 6
            %                 6 3 4];
            %             obj.cases{2}=[1 5 4
            %                 5 2 6
            %                 5 6 4
            %                 6 3 4];
            %             obj.cases{3}=[1 2 5
            %                 1 5 6
            %                 5 3 6
            %                 1 6 4];
            %             obj.cases{4}=[1 2 6
            %                 6 5 4
            %                 2 3 5
            %                 6 2 5];
            %             obj.cases{5}=[1 5 4
            %                 5 6 4
            %                 5 2 6
            %                 2 3 6];
            %             obj.cases{6}=[6 3 4
            %                 5 3 6
            %                 1 5 6
            %                 1 2 5];
            % %             obj.cases{7}=[1 5 8
            % %                 5 2 6
            % %                 5 6 7
            % %                 6 3 7
            % %                 5 7 8
            % %                 8 7 4];
            %             obj.cases{7}=[1 5 8
            %                 5 2 8
            %                 8 2 4
            %                 2 6 4
            %                 4 6 7
            %                 6 3 7];
            
            obj.selectcases = [1 0 0;
                2 0 0;
                3 6 0;
                4 7 0;
                0 5 0;
                0 7 4;
                0 6 3
                0 0 2
                0 0 1];
        end
        
    end
    
    
    
end

