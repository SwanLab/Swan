classdef Filter_P1_LevelSet < Filter_P1 & Filter_LevelSet
    methods
        function obj = Filter_P1_LevelSet(problemID,scale)
            obj@Filter_P1(problemID,scale);
        end
        
        function preProcess(obj)
            preProcess@Filter_P1(obj)
            preProcess@Filter_LevelSet(obj)
        end
        
        function x_gp = getP0fromP1(obj,x)
            if norm(x) == norm(obj.x)
                x_gp=obj.x_reg;
            else
                switch obj.geometry.type
                    case 'TRIANGLE'
                        M2=obj.faireF2(obj.coordinates',obj.connectivities',x);
                    otherwise
                        M2=obj.computeRHS(x);
                end
                x_gp = obj.P_operator*M2;
                obj.x_reg=x_gp;
            end
        end
    end
end

%                                         figure(3)
%                                         hold on
%                                         x_global=coord(connec(cut_elem(ielem),:)',1);
%                                         y_global=coord(connec(cut_elem(ielem),:)',2);
%                                         for i=1:size(x_global,1)
%                                             vec=[i;i+1];
%                                             if i+1 > size(x_global,1)
%                                                 vec(2)=1;
%                                             end
%                                             plot(x_global(vec,:),y_global(vec,:),'b')
%                                         end
%
%                                         for idel=1:size(delaunay_connec,1)
%                                             x_delaunay=mesh_del.coord(delaunay_connec(idel,:)',1);
%                                             y_delaunay=mesh_del.coord(delaunay_connec(idel,:)',2);
%                                             for i=1:size(x_delaunay,1)
%                                                 vec=[i;i+1];
%                                                 if i+1 > size(x_delaunay,1)
%                                                     vec(2)=1;
%                                                 end
%                                                 hold on
%                                                 plot(x_delaunay(vec,:),y_delaunay(vec,:),'-r')
%
%                                                 drawnow
%                                             end
%                                         end
%                                         text(x_global,y_global,num2str(x(connec(cut_elem(ielem),:)')))
%                                         %title(num2str(phi(cut_elem(ielem))))
%                                         close (3)