classdef Quadrilateral_Bilinear < Interpolation

    methods (Access = public)
     function obj = Quadrilateral_Bilinear(mesh)
            obj = obj@Interpolation(mesh);
            obj.type = 'QUADRILATERAL';
            obj.order = 'LINEAR';
            obj.ndime = 2;
            obj.nnode = 4;
            obj.pos_nodes = [-1 -1; 1 -1; 1 1; -1 1];
            obj.dvolu = 4;           
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
        
        function computeShapeDeriv(obj,posgp)
            obj.shape = [];
            obj.deriv = [];
            s = posgp(1,:);
            t = posgp(2,:);
            obj.shape = [(ones(1,size(posgp,2))-t-s+s.*t)*0.25;
                        (ones(1,size(posgp,2))-t+s-s.*t)*0.25;
                        (ones(1,size(posgp,2))+t+s+s.*t)*0.25;
                        (ones(1,size(posgp,2))+t-s-s.*t)*0.25];

            obj.deriv(1,1,:) = (-ones(1,size(posgp,2))+t)*0.25;
            obj.deriv(1,2,:) = (+ones(1,size(posgp,2))-t)*0.25;
            obj.deriv(1,3,:) = (+ones(1,size(posgp,2))+t)*0.25;
            obj.deriv(1,4,:) = (-ones(1,size(posgp,2))-t)*0.25;
            obj.deriv(2,1,:) = (-ones(1,size(posgp,2))+s)*0.25;
            obj.deriv(2,2,:) = (-ones(1,size(posgp,2))-s)*0.25;
            obj.deriv(2,3,:) = (+ones(1,size(posgp,2))+s)*0.25;
            obj.deriv(2,4,:) = (+ones(1,size(posgp,2))-s)*0.25;            
        end

	end
end

