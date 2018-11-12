function [ posgp,weigp ] = quad_gauss_const( ndime,ngaus )
%**** This routine sets up the integration constants of open
%**** integration rules for brick elements:
%
%        NDIME = 1             NDIME = 2             NDIME = 3
%
%    NGAUS  EXACT POL.     NGAUS  EXACT POL.     NGAUS  EXACT POL. 
%    -----  ----------     -----  ----------     -----  ----------
%      1      q1           1 x 1     q1          1x1x1     q1	
%      2      q3           2 x 2     q3          2x2x2     q3   
%      3      q5           3 x 3     q5          3x3x3     q5
%      4      q7           4 x 4     q7          4x4x4     q7
%
%----------------------------------------------------------------------------

if(ndime==1)
   nlocs=ngaus;
elseif(ndime==2)
    switch  ngaus
        case 1
            nlocs = 1;
        case 4
            nlocs = 2;
        case 9
            nlocs = 3;
        otherwise
            error('undifined ngaus')
    end
elseif(ndime==3)
    switch  ngaus
        case 8
            nlocs = 2;
        otherwise
            error('undifined ngaus')
    end
end

if(nlocs==1)
    posgl(1)=0.0;
    weigl(1)=2.0;
elseif(nlocs==2)
    posgl(1)=-0.577350269189626;
    posgl(2)= 0.577350269189626;
    weigl(1)= 1.0;
    weigl(2)= 1.0;
elseif(nlocs==3)
    posgl(1)=-0.774596669241483;
    posgl(2)= 0.0;
    posgl(3)= 0.774596669241483;
    weigl(1)= 0.555555555555556;
    weigl(2)= 0.888888888888889;
    weigl(3)= 0.555555555555556;
elseif(nlocs==4)
    posgl(1)=-0.861136311594053;
    posgl(2)=-0.339981043584856;
    posgl(3)= 0.339981043584856;
    posgl(4)= 0.861136311594053;
    weigl(1)= 0.347854845137454;
    weigl(2)= 0.652145154862546;
	weigl(3)= 0.652145154862546;
	weigl(4)= 0.347854845137454;
end

if(ndime==1)
    igaus=0;
    for ilocs=1:nlocs
        igaus=igaus+1;
        weigp(  igaus)=weigl(ilocs);
        posgp(1,igaus)=posgl(ilocs);
    end
elseif(ndime==2)
    igaus=0;
    for ilocs=1:nlocs
        for jlocs=1:nlocs
            igaus=igaus+1;
            weigp(  igaus)=weigl(ilocs)*weigl(jlocs);
            posgp(1,igaus)=posgl(ilocs);
            posgp(2,igaus)=posgl(jlocs);
        end
    end
elseif(ndime==3)
    igaus=0;
    for ilocs=1:nlocs
        for jlocs=1:nlocs
            for klocs=1:nlocs
                igaus=igaus+1;
                weigp(  igaus)=weigl(ilocs)*weigl(jlocs)*weigl(klocs);
                posgp(1,igaus)=posgl(ilocs);
                posgp(2,igaus)=posgl(jlocs);
                posgp(3,igaus)=posgl(klocs);
            end
        end
    end
end

end

