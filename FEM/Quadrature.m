classdef Quadrature
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        posgp
        weigp
        ngaus
    end
    
    methods
        function obj = Quadrature(geometryType,order)
            switch geometryType
                
                case 'TRIANGLE'
                    
                    switch order
                        case 'CONSTANT'
                            obj.ngaus = 1;          
                            obj.weigp = 1/2;    %revisar
                            obj.posgp = [1/3;1/3];
                        case 'LINEAR' 
                            obj.ngaus = 1;          % Linear triangle
                            obj.weigp = @(igauss) {1/2};
                            obj.posgp = [1/3;1/3];
                        case 'QUADRATIC'
                            obj.ngaus = 3;          % Linear triangle
                            obj.weigp = @(igauss){1/3;1/3;1/3};
                            obj.posgp = [0,0.5;0.5,0;0.5,0.5]';
                        otherwise
                            error('Invalid interpolation order for element TRIANGLE.');
                    end
                case 'Triangle_Linear_Mass'
                    obj.ngaus=3;
                    obj.posgp(1,1)= 2.0/3.0;
                    obj.posgp(2,1)= 1.0/6.0;
                    obj.posgp(1,2)= 1.0/6.0;
                    obj.posgp(2,2)= 2.0/3.0;
                    obj.posgp(1,3)= 1.0/6.0;
                    obj.posgp(2,3)= 1.0/6.0;
                    obj.weigp(  1)= 1.0/6.0;
                    obj.weigp(  2)= 1.0/6.0;
                    obj.weigp(  3)= 1.0/6.0;
                case 'QUAD'
                    switch order
                        case 'LINEAR'
                            obj.ngaus = 4;
                            
                            % Compute WEIGP and POSGP
                            a =  0.577350269189626;
                            obj.posgp(1,:) = [-a,-a];
                            obj.posgp(2,:) = [+a,-a];
                            obj.posgp(3,:) = [-a,+a];
                            obj.posgp(4,:) = [+a,+a];
                            obj.posgp = obj.posgp';
                            obj.weigp = @(igauss) {1,1,1,1};%1*ones(1,obj.ngaus);
                        case 'QUADRATIC' %SERENDIPITY, QUADRILATERAL QUADRATIC NOT IMPLEMENTED
                            obj.ngaus = 9;
                            
                            % Compute WEIGP and POSGP
                            a =  0.77459667;
                            obj.posgp(1,:) = [ 0,+a];
                            obj.posgp(2,:) = [ 0, 0];
                            obj.posgp(3,:) = [+a,+a];
                            obj.posgp(4,:) = [-a,-a];
                            obj.posgp(5,:) = [-a, 0];
                            obj.posgp(6,:) = [+a, 0];
                            obj.posgp(7,:) = [+a,-a];
                            obj.posgp(8,:) = [-a,+a];
                            obj.posgp(9,:) = [ 0,-a];
                            obj.posgp = obj.posgp';
                            obj.weigp = @(igauss) {1,1,1,1};%1*ones(1,obj.ngaus);
                        otherwise
                            error('Invalid interpolation order for element QUADRILATERAL.');
                    end
                case 'TETRAHEDRA'
                    obj.ngaus = 1;          % tetrahedra
                    obj.weigp = 1/6;
                    obj.posgp = [1/4 1/4 1/4];
                case 'HEXAHEDRA'
                    obj.ngaus = 8;
                    %COMPUTE WEIGP AND POSGP
                    
                    nlocs = 2;
                    posgl(1)=-0.577350269189626;
                    posgl(2)= 0.577350269189626;
                    weigl(1)= 1.0;
                    weigl(2)= 1.0;
                    igaus=0;
                    
                    for ilocs=1:nlocs
                        for jlocs=1:nlocs
                            for klocs=1:nlocs
                                igaus=igaus+1;
                                obj.weigp(  igaus)=weigl(ilocs)*weigl(jlocs)*weigl(klocs);
                                obj.posgp(1,igaus)=posgl(ilocs);
                                obj.posgp(2,igaus)=posgl(jlocs);
                                obj.posgp(3,igaus)=posgl(klocs);
                            end
                        end
                    end
                otherwise
                    error('Invalid mesh type.')
            end
        end
        
    end
end
