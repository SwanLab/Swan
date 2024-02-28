classdef Quadrature_Quadrilateral < Quadrature

    methods (Access = public)
        
        function computeQuadrature(obj,order)
            computeQuadrature@Quadrature(obj,order);
            switch order
                case {'CONSTANT','LINEAR'}
                    obj.ngaus = 1;
                    obj.posgp(:,1) = [0,0];
                    obj.weigp = 4;
                    
                case {'QUADRATIC','CUBIC'}
                    obj.ngaus = 4;
                    a =  0.577350269189626;
                    obj.posgp(:,1) = [-a,-a];
                    obj.posgp(:,2) = [+a,-a];
                    obj.posgp(:,3) = [-a,+a];
                    obj.posgp(:,4) = [+a,+a];

                    obj.weigp =  [1,1,1,1];

                case {'ORDER4','ORDER5'}
                    obj.ngaus = 9;
                    weigl = [5/9,8/9,5/9];
                    posgl = [-sqrt(3/5),0,sqrt(3/5)];
                    nlocs = 3;
                    igaus = 0;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            igaus = igaus+1;
                            obj.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
                            obj.posgp(1,igaus) = posgl(ilocs);
                            obj.posgp(2,igaus) = posgl(jlocs);
                        end
                    end

                case {'ORDER6','ORDER7'}
                    obj.ngaus = 16;
                    weigl = [3.47854845137454e-01, 6.52145154862546e-01, 6.52145154862546e-01, 3.47854845137454e-01];
                    posgl = [8.61136311594053e-01, 3.39981043584856e-01, -3.39981043584856e-01, -8.61136311594053e-01];
                    nlocs = 4;
                    igaus = 0;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            igaus = igaus+1;
                            obj.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
                            obj.posgp(1,igaus) = posgl(ilocs);
                            obj.posgp(2,igaus) = posgl(jlocs);
                        end
                    end

                case {'ORDER8','ORDER9'}
                    obj.ngaus = 25;
                    weigl = [2.36926885056189e-01, 4.78628670499366e-01, 5.68888888888889e-01, 4.78628670499366e-01, 2.36926885056189e-01];
                    posgl = [9.06179845938664e-01, 5.38469310105683e-01, 0, -5.38469310105683e-01, -9.06179845938664e-01];
                    nlocs = 5;
                    igaus = 0;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            igaus = igaus+1;
                            obj.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
                            obj.posgp(1,igaus) = posgl(ilocs);
                            obj.posgp(2,igaus) = posgl(jlocs);
                        end
                    end

                case {'ORDER10','ORDER11'}
                    obj.ngaus = 36;
                    weigl = [1.71324492379170e-01, 3.60761573048139e-01, 4.67913934572691e-01, 4.67913934572691e-01, 3.60761573048139e-01, 1.71324492379170e-01];
                    posgl = [9.32469514203152e-01, 6.61209386466265e-01, 2.38619186083197e-01, -2.38619186083197e-01, -6.61209386466265e-01, -9.32469514203152e-01];
                    nlocs = 6;
                    igaus = 0;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            igaus = igaus+1;
                            obj.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
                            obj.posgp(1,igaus) = posgl(ilocs);
                            obj.posgp(2,igaus) = posgl(jlocs);
                        end
                    end

                case {'ORDER12','ORDER13'}
                    obj.ngaus = 49;
                    weigl = [1.29484966168870e-01, 2.79705391489277e-01, 3.81830050505119e-01, 4.17959183673469e-01, 3.81830050505119e-01, 2.79705391489277e-01, 1.29484966168870e-01];
                    posgl = [9.49107912342758e-01, 7.41531185599394e-01, 4.05845151377397e-01, 0, -4.05845151377397e-01, -7.41531185599394e-01, -9.49107912342758e-01];
                    nlocs = 7;
                    igaus = 0;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            igaus = igaus+1;
                            obj.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
                            obj.posgp(1,igaus) = posgl(ilocs);
                            obj.posgp(2,igaus) = posgl(jlocs);
                        end
                    end

                case {'ORDER14','ORDER15'}
                    obj.ngaus = 64;
                    weigl = [1.01228536290377e-01, 2.22381034453374e-01, 3.13706645877887e-01, 3.62683783378362e-01, 3.62683783378362e-01, 3.13706645877887e-01, 2.22381034453374e-01, 1.01228536290377e-01];
                    posgl = [9.60289856497536e-01, 7.96666477413627e-01, 5.25532409916329e-01, 1.83434642495650e-01, -1.83434642495650e-01, -5.25532409916329e-01, -7.96666477413627e-01, -9.60289856497536e-01];
                    nlocs = 8;
                    igaus = 0;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            igaus = igaus+1;
                            obj.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
                            obj.posgp(1,igaus) = posgl(ilocs);
                            obj.posgp(2,igaus) = posgl(jlocs);
                        end
                    end

                case {'ORDER16','ORDER17'}
                    obj.ngaus = 81;
                    weigl = [8.12743883615746e-02, 1.80648160694858e-01, 2.60610696402935e-01, 3.12347077040003e-01, 3.30239355001260e-01, 3.12347077040003e-01, 2.60610696402935e-01, 1.80648160694858e-01, 8.12743883615746e-02];
                    posgl = [9.68160239507626e-01, 8.36031107326636e-01, 6.13371432700590e-01, 3.24253423403809e-01, 00000000000000, -3.24253423403809e-01, -6.13371432700590e-01, -8.36031107326636e-01, -9.68160239507626e-01];
                    nlocs = 9;
                    igaus = 0;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            igaus = igaus+1;
                            obj.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
                            obj.posgp(1,igaus) = posgl(ilocs);
                            obj.posgp(2,igaus) = posgl(jlocs);
                        end
                    end

                case {'ORDER18','ORDER19'}
                    obj.ngaus = 100;
                    weigl = [6.66713443086877e-02, 1.49451349150581e-01, 2.19086362515982e-01, 2.69266719309996e-01, 2.95524224714753e-01, 2.95524224714753e-01, 2.69266719309996e-01, 2.19086362515982e-01, 1.49451349150581e-01, 6.66713443086877e-02];
                    posgl = [9.73906528517172e-01, 8.65063366688985e-01, 6.79409568299024e-01, 4.33395394129247e-01, 1.48874338981631e-01, -1.48874338981631e-01, -4.33395394129247e-01, -6.79409568299024e-01, -8.65063366688985e-01, -9.73906528517172e-01];
                    nlocs = 10;
                    igaus = 0;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            igaus = igaus+1;
                            obj.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
                            obj.posgp(1,igaus) = posgl(ilocs);
                            obj.posgp(2,igaus) = posgl(jlocs);
                        end
                    end

                case {'ORDER20','ORDER21'}
                    obj.ngaus = 121;
                    weigl = [5.56685671161735e-02, 1.25580369464905e-01, 1.86290210927734e-01, 2.33193764591990e-01, 2.62804544510247e-01, 2.72925086777901e-01, 2.62804544510247e-01, 2.33193764591990e-01, 1.86290210927734e-01, 1.25580369464905e-01, 5.56685671161735e-02];
                    posgl = [9.78228658146057e-01, 8.87062599768095e-01, 7.30152005574049e-01, 5.19096129206812e-01, 2.69543155952345e-01, 0, -2.69543155952345e-01, -5.19096129206812e-01, -7.30152005574049e-01, -8.87062599768095e-01, -9.78228658146057e-01];
                    nlocs = 11;
                    igaus = 0;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            igaus = igaus+1;
                            obj.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
                            obj.posgp(1,igaus) = posgl(ilocs);
                            obj.posgp(2,igaus) = posgl(jlocs);
                        end
                    end

                case {'ORDER22','ORDER23'}
                    obj.ngaus = 144;
                    weigl = [4.71753363865118e-02, 1.06939325995318e-01, 1.60078328543346e-01, 2.03167426723066e-01, 2.33492536538355e-01, 2.49147045813403e-01, 2.49147045813403e-01, 2.33492536538355e-01, 2.03167426723066e-01, 1.60078328543346e-01, 1.06939325995318e-01, 4.71753363865118e-02];
                    posgl = [9.81560634246719e-01, 9.04117256370475e-01, 7.69902674194305e-01, 5.87317954286617e-01, 3.67831498998180e-01, 1.25233408511469e-01, -1.25233408511469e-01, -3.67831498998180e-01, -5.87317954286617e-01, -7.69902674194305e-01, -9.04117256370475e-01, -9.81560634246719e-01];
                    nlocs = 12;
                    igaus = 0;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            igaus = igaus+1;
                            obj.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
                            obj.posgp(1,igaus) = posgl(ilocs);
                            obj.posgp(2,igaus) = posgl(jlocs);
                        end
                    end
                        
                otherwise
                    error('Invalid interpolation order for element Quadrilateral.');
            end
        end
    end
end