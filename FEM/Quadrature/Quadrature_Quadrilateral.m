classdef Quadrature_Quadrilateral < Quadrature

    methods (Access = public)
        
        function computeQuadrature(obj,order)
            computeQuadrature@Quadrature(obj,order);
            switch order
                case 'CONSTANT'
                    obj.ngaus = 1;
                    obj.posgp(:,1) = [0,0];
                    obj.weigp = 4;
                    
                case 'LINEAR'
                    obj.ngaus = 4;
                    % Compute WEIGP and POSGP
                    a =  0.577350269189626;
                    obj.posgp(:,1) = [-a,-a];
                    obj.posgp(:,2) = [+a,-a];
                    obj.posgp(:,3) = [-a,+a];
                    obj.posgp(:,4) = [+a,+a];
                    
                    obj.weigp =  [1,1,1,1]';%1*ones(1,obj.ngaus);

                case 'QUADRATIC' %SERENDIPITY, QUADRILATERAL QUADRATIC NOT IMPLEMENTED
                    % Copied from down below
                    obj.ngaus = 9;
                    % Compute WEIGP and POSGP
                    a =  0.77459667;
                    obj.posgp(:,1) = [ 0,+a];
                    obj.posgp(:,2) = [ 0, 0];
                    obj.posgp(:,3) = [+a,+a];
                    obj.posgp(:,4) = [-a,-a];
                    obj.posgp(:,5) = [-a, 0];
                    obj.posgp(:,6) = [+a, 0];
                    obj.posgp(:,7) = [+a,-a];
                    obj.posgp(:,8) = [-a,+a];
                    obj.posgp(:,9) = [ 0,-a];
                    
                    obj.weigp =(4/9)*ones(1,obj.ngaus);

                case 'QUADRATICMSS' %SERENDIPITY, QUADRILATERAL QUADRATIC NOT IMPLEMENTED
                    obj.ngaus = 9;
                    % Compute WEIGP and POSGP
                    a =  0.77459667;
                    obj.posgp(:,1) = [ 0,+a];
                    obj.posgp(:,2) = [ 0, 0];
                    obj.posgp(:,3) = [+a,+a];
                    obj.posgp(:,4) = [-a,-a];
                    obj.posgp(:,5) = [-a, 0];
                    obj.posgp(:,6) = [+a, 0];
                    obj.posgp(:,7) = [+a,-a];
                    obj.posgp(:,8) = [-a,+a];
                    obj.posgp(:,9) = [ 0,-a];
                    
                    obj.weigp =(4/9)*ones(1,obj.ngaus);
                    
                case 'QUADRATICMASS'
                    posgl(1) =-0.774596669241483;
                    posgl(2) = 0.0;
                    posgl(3) = 0.774596669241483;
                    weigl(1) = 0.555555555555556;
                    weigl(2) = 0.888888888888889;
                    weigl(3) = 0.555555555555556;
                    obj.ngaus = 9;
                    igaus = 0;
                    nlocs = 3;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            igaus = igaus+1;
                            obj.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
                            obj.posgp(1,igaus) = posgl(ilocs);
                            obj.posgp(2,igaus) = posgl(jlocs);
                        end
                    end
                    
                case 'ORDER10'
                    obj.ngaus = 24;
                    obj.weigp = [0.09754131164361921;0.2117163509105210;0.2255355359769118;0.09730415294011353;0.2255978016857150;0.3510531956811132;0.3511314245095946;0.2118525411926204;0.2116201646536030;0.3511857026570127;0.3512637749060175;0.2256542961958117;0.09746993474514315;0.2257166005443878;0.2117562971146869;0.09723199711654983;0.06609123516265450;0.06607427278802323;0.06606825552658292;0.06605163800172888;0.04798054519241257;0.04797022666161702;0.04807168904439760;0.04806105514916250];
                    obj.posgp = [-0.1886381798247768,-0.9534228278198672;
                                    0.3158243867065065,-0.8124679583416120;
                                    0.7122535487614264,-0.5253420828029804;
                                    0.9536499381198605,-0.1884848209339622;
                                    -0.5255140441450369,-0.7118387124823607;
                                    -0.04156622116123301,-0.4250108457039333;
                                    0.4249386754351080,-0.04191684210258181;
                                    0.8124112175549880,0.3156521043607226;
                                    -0.8126297846315392,-0.3155908346134177;
                                    -0.4247553230472783,0.04140201954870253;
                                    0.04175248769766685,0.4246831988441449;
                                    0.5251289118559497,0.7121637731013759;
                                    -0.9534285320988584,0.1886914057472521;
                                    -0.7117485896157119,0.5253016177503731;
                                    -0.3154177460532097,0.8125735341734832;
                                    0.1885390832384737,0.9536564551709946;
                                    -0.8257630699887589,-0.9394679906605139;
                                    0.9394356990536500,-0.8259480185291852;
                                    -0.9395368212945202,0.8256048988564928;
                                    0.8257904012838660,0.9395040010720435;
                                    -0.9827093489403464,-0.6980866624366492;
                                    0.6978195696956143,-0.9827223639208844;
                                    -0.6983291127406627,0.9825558709397986;
                                    0.9825693832248631,0.6980615873134067]';
                        
                otherwise
                    error('Invalid interpolation order for element Quadrilateral.');
            end
        end
    end
end