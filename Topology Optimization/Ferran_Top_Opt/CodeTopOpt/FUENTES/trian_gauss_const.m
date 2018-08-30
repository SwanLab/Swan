function [ posgp,weigp ] = trian_gauss_const( ndime,ngaus )

%
%**** This routine sets up the integration constants of open rules for
%**** triangles and tetrahedra
%
%            NDIME = 2             NDIME = 3
%
%         NGAUS  EXACT POL.     NGAUS  EXACT POL. 
%         -----  ----------     -----  ----------
%           1       p1            1       p1
%           3       p2            4       p2
%           4       p3            5       p3
%           6       p4           11       p4
%           7       p5           14       p5
%          13       p9
%
%***********************************************************************

%
%*** Line integral (the same as for brick elements)
%
istop=0;
if(ndime==1)
    if(ngaus==1)
        posgp(1,1)=0.0;
        weigp(  1)=2.0;
    elseif(ngaus==2)
        posgp(1,1)=-0.577350269189626;
        posgp(1,2)= 0.577350269189626;
        weigp(  1)= 1.0;
        weigp(  2)= 1.0;
    elseif(ngaus==3)
        posgp(1,1)=-0.774596669241483;
        posgp(1,2)= 0.0;
        posgp(1,3)= 0.774596669241483;
        weigp(  1)= 0.555555555555556;
        weigp(  2)= 0.888888888888889;
        weigp(  3)= 0.555555555555556;
    else
        istop=1;
    end
end
        
%        
%*** Area integral (triangles)
%
if(ndime==2)
    if(ngaus==1)
        posgp(1,1)= 1.0/3.0;
        posgp(2,1)= 1.0/3.0;
        weigp(  1)= 1.0/2.0;
    elseif(ngaus==3)
        posgp(1,1)= 2.0/3.0;
        posgp(2,1)= 1.0/6.0;
        posgp(1,2)= 1.0/6.0;
        posgp(2,2)= 2.0/3.0;
        posgp(1,3)= 1.0/6.0;
        posgp(2,3)= 1.0/6.0;
        weigp(  1)= 1.0/6.0;
        weigp(  2)= 1.0/6.0;
        weigp(  3)= 1.0/6.0;
    elseif(ngaus==4)
        posgp(1,1)= 1.0/3.0;
        posgp(2,1)= 1.0/3.0;
        posgp(1,2)= 1.0/5.0;
        posgp(2,2)= 1.0/5.0;
        posgp(1,3)= 3.0/5.0;
        posgp(2,3)= 1.0/5.0;
        posgp(1,4)= 1.0/5.0;
        posgp(2,4)= 3.0/5.0;
        weigp(  1)=-27.0/96.0;
        weigp(  2)= 25.0/96.0;
        weigp(  3)= 25.0/96.0;
        weigp(  4)= 25.0/96.0;
    elseif(ngaus==6)
        ex1 = 0.816847572980459;
        et1 = 0.091576213509771;
        ez1 = 0.091576213509771;
        ex2 = 0.108103018168070;
        et2 = 0.445948490915965;
        ez2 = 0.445948490915965;
        posgp(1,1)= ex1;
        posgp(2,1)= et1;
        posgp(1,2)= et1;
        posgp(2,2)= ez1;
        posgp(1,3)= ez1;
        posgp(2,3)= ex1;
        posgp(1,4)= ex2;
        posgp(2,4)= et2;
        posgp(1,5)= et2;
        posgp(2,5)= ez2;
        posgp(1,6)= ez2;
        posgp(2,6)= ex2;
        a = 0.054975870996713638;
        b = 0.1116907969117165;
        weigp(1)  = a;
        weigp(2)  = a;
        weigp(3)  = a;
        weigp(4)  = b;
        weigp(5)  = b;
        weigp(6)  = b;
    elseif(ngaus==7)
        posgp(1,1)= 0.333333333333333;
        posgp(2,1)= 0.333333333333333;
        posgp(1,2)= 0.059715871789770;
        posgp(2,2)= 0.470142064105115;
        posgp(1,3)= 0.470142064105115;
        posgp(2,3)= 0.059715871789770;
        posgp(1,4)= 0.470142064105115;
        posgp(2,4)= 0.470142064105115;
        posgp(1,5)= 0.797426985353087;
        posgp(2,5)= 0.101286507323456;
        posgp(1,6)= 0.101286507323456;
        posgp(2,6)= 0.797426985353087;
        posgp(1,7)= 0.101286507323456;
        posgp(2,7)= 0.101286507323456;
        weigp(  1)= 0.112500000000000;
        weigp(  2)= 0.066197074949741;
        weigp(  3)= 0.066197074949741;
        weigp(  4)= 0.066197074949741;
        weigp(  5)= 0.062969587743282;
        weigp(  6)= 0.062969587743282;
        weigp(  7)= 0.062969587743282;
    elseif(ngaus==13)
        a = 0.333333333333333;
        b = 0.479308067841920;
        c = 0.869739794195568;
        d = 0.638444188569810;
        e = 0.260345966079040;
        f = 0.065130102902216;
        g = 0.312865496004874;
        h = 0.048690315425316;
        w1=-0.149570044467670/2.0;
        w2= 0.175615257433204/2.0;
        w3= 0.053347235608839/2.0;
        w4= 0.077113760890257/2.0;
        posgp(1, 1)= a;
        posgp(2, 1)= a;
        posgp(1, 2)= e;
        posgp(2, 2)= e;
        posgp(1, 3)= b;
        posgp(2, 3)= e;
        posgp(1, 4)= e;
        posgp(2, 4)= b;
        posgp(1, 5)= f;
        posgp(2, 5)= f;
        posgp(1, 6)= c;
        posgp(2, 6)= f;
        posgp(1, 7)= f;
        posgp(2, 7)= c;
        posgp(1, 8)= d;
        posgp(2, 8)= g;
        posgp(1, 9)= d;
        posgp(2, 9)= h;
        posgp(1,10)= g;
        posgp(2,10)= d;
        posgp(1,11)= g;
        posgp(2,11)= h;
        posgp(1,12)= h;
        posgp(2,12)= d;
        posgp(1,13)= h;
        posgp(2,13)= g;
        weigp( 1) = w1;
        weigp( 2) = w2;
        weigp( 3) = w2;
        weigp( 4) = w2;
        weigp( 5) = w3;
        weigp( 6) = w3;
        weigp( 7) = w3;
        weigp( 8) = w4;
        weigp( 9) = w4;
        weigp(10) = w4;
        weigp(11) = w4;
        weigp(12) = w4;
        weigp(13) = w4;
    else
        istop=1;
    end
end
%
%*** Volume integral ( tetrahedra )
%
if(ndime==3)
    if(ngaus==1)
        posgp(1,1)= 1.0/4.0;
        posgp(2,1)= 1.0/4.0;
        posgp(3,1)= 1.0/4.0;
        weigp(1)  = 1.0/6.0;
    elseif(ngaus==4)
        a=0.5854101966249685;
        b=0.1381966011250105;
        posgp(1,1)= b;
        posgp(2,1)= b;
        posgp(3,1)= b;
        posgp(1,2)= a;
        posgp(2,2)= b;
        posgp(3,2)= b;
        posgp(1,3)= b;
        posgp(2,3)= a;
        posgp(3,3)= b;
        posgp(1,4)= b;
        posgp(2,4)= b;
        posgp(3,4)= a;
        weigp(  1)= 1.0/24.0;
        weigp(  2)= 1.0/24.0;
        weigp(  3)= 1.0/24.0;
        weigp(  4)= 1.0/24.0;
    elseif(ngaus==5)
        posgp(1,1)= 1.0/4.0;
        posgp(2,1)= 1.0/4.0;
        posgp(3,1)= 1.0/4.0;
        posgp(1,2)= 1.0/6.0;
        posgp(2,2)= 1.0/6.0;
        posgp(3,2)= 1.0/6.0;
        posgp(1,3)= 1.0/2.0;
        posgp(2,3)= 1.0/6.0;
        posgp(3,3)= 1.0/6.0;
        posgp(1,4)= 1.0/6.0;
        posgp(2,4)= 1.0/2.0;
        posgp(3,4)= 1.0/6.0;
        posgp(1,5)= 1.0/6.0;
        posgp(2,5)= 1.0/6.0;
        posgp(3,5)= 1.0/2.0;
        weigp(  1)=-2.0/15.0;
        weigp(  2)= 1.5/20.0;
        weigp(  3)= 1.5/20.0;
        weigp(  4)= 1.5/20.0;
        weigp(  5)= 1.5/20.0;
    elseif(ngaus==11)
        a=0.3994035761667992;
        b=0.1005964238332008;
        c=343.0/7500.0/6.0;
        d=56.0/375.0/6.0;
        posgp(1,1) = 1.0/4.0;
        posgp(2,1) = 1.0/4.0;
        posgp(3,1) = 1.0/4.0;
        posgp(1,2) = 11.0/14.0;
        posgp(2,2) = 1.0/14.0;
        posgp(3,2) = 1.0/14.0;
        posgp(1,3) = 1.0/14.0;
        posgp(2,3) = 11.0/14.0;
        posgp(3,3) = 1.0/14.0;
        posgp(1,4) = 1.0/14.0;
        posgp(2,4) = 1.0/14.0;
        posgp(3,4) = 11.0/14.0;
        posgp(1,5) = 1.0/14.0;
        posgp(2,5) = 1.0/14.0;
        posgp(3,5) = 1.0/14.0;
        posgp(1,6) = a;
        posgp(2,6) = a;
        posgp(3,6) = b;
        posgp(1,7) = a;
        posgp(2,7) = b;
        posgp(3,7) = a;
        posgp(1,8) = a;
        posgp(2,8) = b;
        posgp(3,8) = b;
        posgp(1,9) = b;
        posgp(2,9) = a;
        posgp(3,9) = a;
        posgp(1,10)= b;
        posgp(2,10)= a;
        posgp(3,10)= b;
        posgp(1,11)= b;
        posgp(2,11)= b;
        posgp(3,11)= a;
        weigp(1)   =-148.0/1875.0/6.0;
        weigp(2)   = c;
        weigp(3)   = c;
        weigp(4)   = c;
        weigp(5)   = c;
        weigp(6)   = d;
        weigp(7)   = d;
        weigp(8)   = d;
        weigp(9)   = d;
        weigp(10)  = d;
        weigp(11)  = d;
    elseif(ngaus==14)
        a=0.0673422422100983;
        b=0.3108859192633005;
        c=0.7217942490673264;
        d=0.0927352503108912;
        e=0.4544962958743506;
        f=0.0455037041256494;
        p=0.1126879257180162/6.0;
        q=0.0734930431163619/6.0;
        r=0.0425460207770812/6.0;
        posgp(1,1) = a;
        posgp(2,1) = b;
        posgp(3,1) = b;
        posgp(1,2) = b;
        posgp(2,2) = a;
        posgp(3,2) = b;
        posgp(1,3) = b;
        posgp(2,3) = b;
        posgp(3,3) = a;
        posgp(1,4) = b;
        posgp(2,4) = b;
        posgp(3,4) = b;
        posgp(1,5) = c;
        posgp(2,5) = d;
        posgp(3,5) = d;
        posgp(1,6) = d;
        posgp(2,6) = c;
        posgp(3,6) = d;
        posgp(1,7) = d;
        posgp(2,7) = d;
        posgp(3,7) = c;
        posgp(1,8) = d;
        posgp(2,8) = d;
        posgp(3,8) = d;
        posgp(1,9) = e;
        posgp(2,9) = e;
        posgp(3,9) = f;
        posgp(1,10)= e;
        posgp(2,10)= f;
        posgp(3,10)= e;
        posgp(1,11)= e;
        posgp(2,11)= f;
        posgp(3,11)= f;
        posgp(1,12)= f;
        posgp(2,12)= e;
        posgp(3,12)= e;
        posgp(1,13)= f;
        posgp(2,13)= e;
        posgp(3,13)= f;
        posgp(1,14)= f;
        posgp(2,14)= f;
        posgp(3,14)= e;
        weigp(1)   = p;
        weigp(2)   = p;
        weigp(3)   = p;
        weigp(4)   = p;
        weigp(5)   = q;
        weigp(6)   = q;
        weigp(7)   = q;
        weigp(8)   = q;
        weigp(9)   = r;
        weigp(10)  = r;
        weigp(11)  = r;
        weigp(12)  = r;
        weigp(13)  = r;
        weigp(14)  = r;
    else
        istop=1;
    end
end
%      
%*** Errors
%
if(istop==1)
    fprintf (1,'TRIANG_GAUSS_CONST:  NOT AVAILABLE QUADRATURE\n');
end

end

