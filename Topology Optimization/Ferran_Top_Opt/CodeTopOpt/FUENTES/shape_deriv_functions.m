function [shape,deriv,heslo] = shape_deriv_functions(igaus,posgp,ptype,etype,nnode,neres)
%**** This routine evaluates shape functions and their first and
%**** second derivatives for 2-d continuous standar interpolation 
%**** elements.
%
%     TRIANGLES       3   4   6  7  &  10  nodes
%     QUADRILATERALS  4   5   8  9  &  16  nodes
%     HEXAHEDRA       8   
%
%     s : xi coordinate 
%     t : eta coordinate
%     u : zeta coordinate (for 3D)
%     element.type : element type (quadrilaterals, triangles);
%     nnode : number of nodes by element
%     neres : if 1, computes hessian matrix
%***********************************************************************

switch ptype
      
    case '2D'
        ndime = 2; 
        s=posgp(1,igaus);t=posgp(2,igaus);  
        deriv = zeros(ndime,nnode);
        shape = zeros(1,nnode);
        heslo = zeros(3,nnode);
        switch  etype
            case 'TRIANGLE'
                if(nnode==3)
                    shape(1)=1.0-s-t;
                    shape(2)=s;
                    shape(3)=t;                                          %  3
                    deriv(1,1)=-1.0;                                     %
                    deriv(1,2)= 1.0;                                     %
                    deriv(1,3)= 0.0;                                     %
                    deriv(2,1)=-1.0;                                     %  1       2
                    deriv(2,2)= 0.0;
                    deriv(2,3)= 1.0;
                elseif(nnode==4)
                    a1=1.0-s-t;
                    a2=s;
                    a3=t;                                                %  3
                    shape(1)=a1;                                         %
                    shape(2)=a2;                                         %
                    shape(3)=a3;                                         %     4
                    shape(4)=27.0*a1*a2*a3;                              %
                    deriv(1,1)=-1.0;                                     %  1        2
                    deriv(1,2)= 1.0;
                    deriv(1,3)= 0.0;
                    deriv(1,4)= 27.0*a3*(a1-a2);
                    deriv(2,1)=-1.0;
                    deriv(2,2)= 0.0;
                    deriv(2,3)= 1.0;
                    deriv(2,4)= 27.0*a2*(a1-a3);
                    if(neres==1)
                        heslo(1,4)=-54.0*a3;
                        heslo(2,4)=-54.0*a2;
                        heslo(3,4)= 54.0*a1-27.0;
                    end
                    for i=1:3
                        shape(  i)=shape(  i)-1./3.*shape(  4);
                        deriv(1,i)=deriv(1,i)-1./3.*deriv(1,4);
                        deriv(2,i)=deriv(2,i)-1./3.*deriv(2,4);
                        if (neres==1);
                            heslo(1,i)=heslo(1,i)-1./3.*heslo(1,4);
                            heslo(2,i)=heslo(2,i)-1./3.*heslo(2,4);
                            heslo(3,i)=heslo(3,i)-1./3.*heslo(3,4);
                        end
                    end
                elseif(nnode==6)
                    a1=1.0-s-t;
                    a2=s;
                    a3=t;
                    shape( 1)=(2.0*a1-1.0)*a1;   %  3
                    shape( 2)=(2.0*a2-1.0)*a2;   %
                    shape( 3)=(2.0*a3-1.0)*a3;   %
                    shape( 4)=4.0*a1*a2;         %  6      5
                    shape( 5)=4.0*a2*a3;         %
                    shape( 6)=4.0*a1*a3;         %
                    deriv(1,1)= 1.0-4.0*a1;      %  1     4     2
                    deriv(1,2)= 4.0*a2-1.0;
                    deriv(1,3)= 0.0;
                    deriv(1,4)= 4.0*(a1-a2);
                    deriv(1,5)= 4.0*a3;
                    deriv(1,6)=-4.0*a3;
                    deriv(2,1)= 1.0-4.0*a1;
                    deriv(2,2)= 0.0;
                    deriv(2,3)= 4.0*a3-1.0;
                    deriv(2,4)=-4.0*a2;
                    deriv(2,5)= 4.0*a2;
                    deriv(2,6)= 4.0*(a1-a3);
                    if(neres==1)
                        heslo(1,1)= 4.0;
                        heslo(1,2)= 4.0;
                        heslo(1,4)=-8.0;
                        heslo(2,1)= 4.0;
                        heslo(2,3)= 4.0;
                        heslo(2,6)=-8.0;
                        heslo(3,1)= 4.0;
                        heslo(3,4)=-4.0;
                        heslo(3,5)= 4.0;
                        heslo(3,6)=-4.0;
                    end
                elseif(nnode==7)
                    a1=1.0-s-t;
                    a2=s;
                    a3=t;
                    v123=-1.0/9.0;
                    v456= 4.0/9.0;                        %
                    shape( 1)=(2.0*a1-1.0)*a1;            %  3
                    shape( 2)=(2.0*a2-1.0)*a2;            %
                    shape( 3)=(2.0*a3-1.0)*a3;            %
                    shape( 4)=4.0*a1*a2;                  %  6      5
                    shape( 5)=4.0*a2*a3;                  %      7
                    shape( 6)=4.0*a1*a3;                  %  1      4     2
                    shape( 7)=27.0*a1*a2*a3;              % bubble function
                    for inode=1:3
                        shape(inode  )=shape(inode  )-v123*shape(7);
                        shape(inode+3)=shape(inode+3)-v456*shape(7);
                    end
                    deriv(1, 1)= 1.0-4.0*a1;
                    deriv(1, 2)= 4.0*a2-1.0;
                    deriv(1, 3)= 0.0;
                    deriv(1, 4)= 4.0*(a1-a2);
                    deriv(1, 5)= 4.0*a3;
                    deriv(1, 6)=-4.0*a3;
                    deriv(1, 7)= 27.0*a3*(a1-a2);
                    deriv(2, 1)= 1.0-4.0*a1;
                    deriv(2, 2)= 0.0;
                    deriv(2, 3)= 4.0*a3-1.0;
                    deriv(2, 4)=-4.0*a2;
                    deriv(2, 5)= 4.0*a2;
                    deriv(2, 6)= 4.0*(a1-a3);
                    deriv(2, 7)= 27.0*a2*(a1-a3);
                    for inode=1:3
                        deriv(1,inode  )=deriv(1,inode  )-v123*deriv(1,7);
                        deriv(1,inode+3)=deriv(1,inode+3)-v456*deriv(1,7);
                        deriv(2,inode  )=deriv(2,inode  )-v123*deriv(2,7);
                        deriv(2,inode+3)=deriv(2,inode+3)-v456*deriv(2,7);
                    end
                    if(neres==1)
                        heslo(1,1)= 4.0;
                        heslo(1,2)= 4.0;
                        heslo(1,4)=-8.0;
                        heslo(1,7)=-54.0*a3;
                        heslo(2,1)= 4.0;
                        heslo(2,3)= 4.0;
                        heslo(2,6)=-8.0;
                        heslo(2,7)=-54.0*a2;
                        heslo(3,1)= 4.0;
                        heslo(3,4)=-4.0;
                        heslo(3,5)= 4.0;
                        heslo(3,6)=-4.0;
                        heslo(3,7)= 54.0*a1-27.0;
                        for inode=1:3
                            heslo(1,inode  )=heslo(1,inode  )-v123*heslo(1,7);
                            heslo(1,inode+3)=heslo(1,inode+3)-v456*heslo(1,7);
                            heslo(2,inode  )=heslo(2,inode  )-v123*heslo(2,7);
                            heslo(2,inode+3)=heslo(2,inode+3)-v456*heslo(2,7);
                            heslo(3,inode  )=heslo(3,inode  )-v123*heslo(3,7);
                            heslo(3,inode+3)=heslo(3,inode+3)-v456*heslo(3,7);
                        end
                    end
                elseif(nnode==10)
                    c=9.0/2.0;
                    a1=1.0-s-t;
                    a2=2.0/3.0-s-t;
                    a3=1.0/3.0-s-t;
                    shape( 1)=c*a1*a2*a3;                                %  3
                    shape( 2)=c*(1.0/3.0-s)*(2.0/3.0-s)*s;               %
                    shape( 3)=c*(1.0/3.0-t)*(2.0/3.0-t)*t;               %
                    shape( 4)= 3.0*c*a1*a2*s;                            %  8    7
                    shape( 5)=-3.0*c*a1*(1.0/3.0-s);                     %
                    shape( 6)=-3.0*c*(1.0/3.0-s)*s*t;                    %
                    shape( 7)=-3.0*c*s*(1.0/3.0-t)*t;                    %  9   10    6
                    shape( 8)=-3.0*c*a1*(1.0/3.0-t)*t;                   %
                    shape( 9)= 3.0*c*a1*a2*t;                            %
                    shape(10)= 6.0*c*a1*s*t;                             %  1    4    5    2
                    deriv(1, 1)=-c*(a1*a2+a1*a3+a2*a3);
                    deriv(1, 2)=-c*((2.0/3.0-s)*s+(1.0/3.0-s)*s-(1.0/3.0-s)*(2.0/3.0-s));
                    deriv(1, 3)=0.0;
                    deriv(1, 4)= 3.0*c*(a1*a2-a1*s-a2*s);
                    deriv(1, 5)=-3.0*c*(a1*(1.0/3.0-s)-a1*s-(1.0/3.0-s)*s);
                    deriv(1, 6)=-3.0*c*((1.0/3.0-s)*t-s*t);
                    deriv(1, 7)=-3.0*c*((1.0/3.0-t)*t);
                    deriv(1, 8)= 3.0*c*((1.0/3.0-t)*t);
                    deriv(1, 9)= 3.0*c*(-a1*t-a2*t);
                    deriv(1,10)= 6.0*c*(a1*t-s*t);
                    deriv(2, 1)=-c*(a1*a2+a1*a3+a2*a3);
                    deriv(2, 2)= 0.0;
                    deriv(2, 3)=-c*((2.0/3.0-t)*t+(1.0/3.0-t)*t-(1.0/3.0-t)*(2.0/3.0-t));
                    deriv(2, 4)= 3.0*c*(-a1*s-a2*s);
                    deriv(2, 5)=-3.0*c*(-(1.0/3.0-s)*s);
                    deriv(2, 6)=-3.0*c*((1.0/3.0-s)*s);
                    deriv(2, 7)=-3.0*c*((1.0/3.0-t)*s-s*t);
                    deriv(2, 8)=-3.0*c*(-(1.0/3.0-t)*t-a1*t+a1*(1.0/3.0-t));
                    deriv(2, 9)= 3.0*c*(-a1*t-a2*t+a1*a2);
                    deriv(2,10)= 6.0*c*(a1*s-s*t);
                    if (neres==1)
                        heslo(1, 1)= 2.0*c*(a1+a2+a3);
                        heslo(1, 2)=-c*(2.0-6.0*s);
                        heslo(1, 3)= 0.0;
                        heslo(1, 4)= 3.0*c*(-2.0*a1-2.0*a2+2.0*s);
                        heslo(1, 5)=-3.0*c*(-2.0*a1-2.0*(1.0/3.0-s)+2.0*s);
                        heslo(1, 6)= 3.0*c*2.0*t;
                        heslo(1, 7)= 0.0;
                        heslo(1, 8)= 0.0;
                        heslo(1, 9)= 3.0*c*2.0*t;
                        heslo(1,10)=-6.0*c*2.0*t;
                        heslo(2, 1)= c*(2.0*a1+2.0*a2+2.0*a3);
                        heslo(2, 2)= 0.0;
                        heslo(2, 3)=-c*(2.0-6.0*t);
                        heslo(2, 4)= 3.0*c*2.0*s;
                        heslo(2, 5)= 0.0;
                        heslo(2, 6)= 0.0;
                        heslo(2, 7)= 3.0*c*2.0*s;
                        heslo(2, 8)=-3.0*c*(-2.0*a1-2.0*(1.0/3.0-t)+2.0*t);
                        heslo(2, 9)= 3.0*c*(-2.0*a1-2.0*a2+2.0*t);
                        heslo(2,10)=-6.0*c*2.0*s;
                        heslo(3, 1)= 2.0*c*(a1+a2+a3);
                        heslo(3, 2)= 0.0;
                        heslo(3, 3)= 0.0;
                        heslo(3, 4)= 3.0*c*(-a1-a2+2.0*s);
                        heslo(3, 5)=-3.0*c*(-(1.0/3.0-s)+s);
                        heslo(3, 6)=-3.0*c*(1.0/3.0-2.0*s);
                        heslo(3, 7)=-3.0*c*(1.0/3.0-2.0*t);
                        heslo(3, 8)= 3.0*c*(1.0/3.0-2.0*t);
                        heslo(3, 9)= 3.0*c*(-a1-a2+2.0*t);
                        heslo(3,10)= 6.0*c*(a1-s-t);
                    end
                end
            case 'QUAD'
                if(nnode==4)
                    st=s*t;
                    shape(1)=(1.-t-s+st)*0.25;           %  4         3
                    shape(2)=(1.-t+s-st)*0.25;           %
                    shape(3)=(1.+t+s+st)*0.25;           %
                    shape(4)=(1.+t-s-st)*0.25;           %
                    deriv(1,1)=(-1.+t)*0.25;             %  1         2
                    deriv(1,2)=(+1.-t)*0.25;
                    deriv(1,3)=(+1.+t)*0.25;
                    deriv(1,4)=(-1.-t)*0.25;
                    deriv(2,1)=(-1.+s)*0.25;
                    deriv(2,2)=(-1.-s)*0.25;
                    deriv(2,3)=(+1.+s)*0.25;
                    deriv(2,4)=(+1.-s)*0.25;
                    if(neres==1)
                        heslo(3,1)= 0.25;
                        heslo(3,2)=-0.25;
                        heslo(3,3)= 0.25;
                        heslo(3,4)=-0.25;
                    end
                elseif(nnode==5)
                    ss=s*s;
                    tt=t*t;
                    st=s*t;
                    t2=2.0*t;
                    s2=2.0*s;
                    shape(1)=(1.-t-s+st)*0.25;        %  4         3
                    shape(2)=(1.-t+s-st)*0.25;        %
                    shape(3)=(1.+t+s+st)*0.25;        %       5
                    shape(4)=(1.+t-s-st)*0.25;        %
                    shape(5)=(1.0-ss)*(1.0-tt);       %  1         2
                    deriv(1,1)=(-1.+t)*0.25;
                    deriv(1,2)=(+1.-t)*0.25;
                    deriv(1,3)=(+1.+t)*0.25;
                    deriv(1,4)=(-1.-t)*0.25;
                    deriv(1,5)=-s2*(1.0-tt);
                    deriv(2,1)=(-1.+s)*0.25;
                    deriv(2,2)=(-1.-s)*0.25;
                    deriv(2,3)=(+1.+s)*0.25;
                    deriv(2,4)=(+1.-s)*0.25;
                    deriv(2,5)=-t2*(1.0-ss);
                    if(neres==1)
                        heslo(1,5)=-2.0*(1.0-tt);
                        heslo(2,5)=-2.0*(1.0-ss);
                        heslo(3,1)= 0.25;
                        heslo(3,2)=-0.25;
                        heslo(3,3)= 0.25;
                        heslo(3,4)=-0.25;
                        heslo(3,5)= s2*t2;
                    end
                    for i=1:4
                        shape(  i)=shape(  i)-0.25*shape(  5);
                        deriv(1,i)=deriv(1,i)-0.25*deriv(1,5);
                        deriv(2,i)=deriv(2,i)-0.25*deriv(2,5);
                        if (neres==1); then
                            heslo(1,i)=heslo(1,i)-0.25*heslo(1,5);
                            heslo(2,i)=heslo(2,i)-0.25*heslo(2,5);
                            heslo(3,i)=heslo(3,i)-0.25*heslo(3,5);
                        end
                    end
                    
                elseif(nnode==8)
                    s2=s*2.0;
                    t2=t*2.0;
                    ss=s*s;
                    tt=t*t;
                    st=s*t;
                    sst=s*s*t;
                    stt=s*t*t;
                    st2=s*t*2.0;
                    shape(1)=(-1.0+st+ss+tt-sst-stt)/4.0;
                    shape(2)=(-1.0-st+ss+tt-sst+stt)/4.0;
                    shape(3)=(-1.0+st+ss+tt+sst+stt)/4.0;
                    shape(4)=(-1.0-st+ss+tt+sst-stt)/4.0;
                    shape(5)=(1.0-t-ss+sst)/2.0;                         %  4    7    3
                    shape(6)=(1.0+s-tt-stt)/2.0;                         %
                    shape(7)=(1.0+t-ss-sst)/2.0;                         %  8         6
                    shape(8)=(1.0-s-tt+stt)/2.0;                         %
                    deriv(1,1)=(t+s2-st2-tt)/4.0;                        %  1    5    2
                    deriv(1,2)=(-t+s2-st2+tt)/4.0;
                    deriv(1,3)=(t+s2+st2+tt)/4.0;
                    deriv(1,4)=(-t+s2+st2-tt)/4.0;
                    deriv(1,5)=-s+st;
                    deriv(1,6)=(1.0-tt)/2.0;
                    deriv(1,7)=-s-st;
                    deriv(1,8)=(-1.0+tt)/2.0;
                    deriv(2,1)=(s+t2-ss-st2)/4.0;
                    deriv(2,2)=(-s+t2-ss+st2)/4.0;
                    deriv(2,3)=(s+t2+ss+st2)/4.0;
                    deriv(2,4)=(-s+t2+ss-st2)/4.0;
                    deriv(2,5)=(-1.0+ss)/2.0;
                    deriv(2,6)=-t-st;
                    deriv(2,7)=(1.0-ss)/2.0;
                    deriv(2,8)=-t+st;
                    if(neres==1)
                        heslo(1,1)= (2.0-t2)/4.0;
                        heslo(1,2)= (2.0-t2)/4.0;
                        heslo(1,3)= (2.0+t2)/4.0;
                        heslo(1,4)= (2.0+t2)/4.0;
                        heslo(1,5)=-1.0+t;
                        heslo(1,7)=-1.0-t;
                        heslo(2,1)= (2.0-s2)/4.0;
                        heslo(2,2)= (2.0+s2)/4.0;
                        heslo(2,3)= (2.0+s2)/4.0;
                        heslo(2,4)= (2.0-s2)/4.0;
                        heslo(2,6)=-1.0-s;
                        heslo(2,8)=-1.0+s;
                        heslo(3,1)= (1.0-s2-t2)/4.0;
                        heslo(3,2)= (-1.0-s2+t2)/4.0;
                        heslo(3,3)= (1.d00+s2+t2)/4.0;
                        heslo(3,4)= (-1.0+s2-t2)/4.0;
                        heslo(3,5)= s;
                        heslo(3,6)=-t;
                        heslo(3,7)=-s;
                        heslo(3,8)= t;
                    end
                elseif(nnode==9)
                    ss=s*s;
                    st=s*t;
                    tt=t*t;
                    s1=s+1.0;
                    t1=t+1.0;
                    s2=s*2.0;
                    t2=t*2.0;
                    s9=s-1.0;
                    t9=t-1.0;                                            %  4      7      3
                    shape( 1)=0.25*s9*st*t9;                             %
                    shape( 2)=0.25*s1*st*t9;                             %
                    shape( 3)=0.25*s1*st*t1;                             %
                    shape( 4)=0.25*s9*st*t1;                             %  8      9      6
                    shape( 5)=0.5*(1.0-ss)*t*t9;                         %
                    shape( 6)=0.5*s*s1*(1.0-tt);                         %
                    shape( 7)=0.5*(1.0-ss)*t*t1;                         %
                    shape( 8)=0.5*s*s9*(1.0-tt);                         %  1      5      2
                    shape( 9)=(1.0-ss)*(1.0-tt);
                    deriv(1,1)= 0.25*t*t9*(-1.0+s2);
                    deriv(1,2)= 0.25*(1.0+s2)*t*t9;
                    deriv(1,3)= 0.25*(1.0+s2)*t*t1;
                    deriv(1,4)= 0.25*(-1.0+s2)*t*t1;
                    deriv(1,5)=-st*t9;
                    deriv(1,6)= 0.5*(1.0+s2)*(1.0-tt);
                    deriv(1,7)=-st*t1;
                    deriv(1,8)= 0.5*(-1.0+s2)*(1.0-tt);
                    deriv(1,9)=-s2*(1.0-tt);
                    deriv(2,1)= 0.25*(-1.0+t2)*s*s9;
                    deriv(2,2)= 0.25*s*s1*(-1.0+t2);
                    deriv(2,3)= 0.25*s*s1*(1.0+t2);
                    deriv(2,4)= 0.25*s*s9*(1.0+t2);
                    deriv(2,5)= 0.5*(1.0-ss)*(-1.0+t2);
                    deriv(2,6)=-st*s1;
                    deriv(2,7)= 0.5*(1.0-ss)*(1.0+t2);
                    deriv(2,8)=-st*s9;
                    deriv(2,9)=-t2*(1.0-ss);
                    if(neres==1)
                        heslo(1,1)= 0.5*t*t9;
                        heslo(1,2)= 0.5*t*t9;
                        heslo(1,3)= 0.5*t*t1;
                        heslo(1,4)= 0.5*t*t1;
                        heslo(1,5)=-t*t9;
                        heslo(1,6)= 1.0-tt;
                        heslo(1,7)=-t*t1;
                        heslo(1,8)= 1.0-tt;
                        heslo(1,9)=-2.0*(1.0-tt);
                        heslo(2,1)= 0.5*s*s9;
                        heslo(2,2)= 0.5*s*s1;
                        heslo(2,3)= 0.5*s*s1;
                        heslo(2,4)= 0.5*s*s9;
                        heslo(2,5)= 1.0-ss;
                        heslo(2,6)=-s*s1;
                        heslo(2,7)= 1.0-ss;
                        heslo(2,8)=-s*s9;
                        heslo(2,9)=-2.0*(1.0-ss);
                        heslo(3,1)= 0.25*(-1.0+t2)*(s9+s);
                        heslo(3,2)= 0.25*(-1.0+t2)*(s1+s);
                        heslo(3,3)= 0.25*(1.0+t2)*(s1+s);
                        heslo(3,4)= 0.25*(1.0+t2)*(s9+s);
                        heslo(3,5)=-s*(-1.0+t2);
                        heslo(3,6)=-t*s1-st;
                        heslo(3,7)=-s*(1.0+t2);
                        heslo(3,8)=-t*s9-st;
                        heslo(3,9)= s2*t2;
                    end
                    
                elseif(nnode==16)
                    a =81.0/256.0;
                    c =1.0/3.0;
                    s1=1.0+s;
                    s2=c+s;
                    s3=c-s;
                    s4=1.0-s;
                    t1=1.0+t;
                    t2=c+t;
                    t3=c-t;
                    t4=1.0-t;
                    shape( 1)=   a*s2*s3*s4*t2*t3*t4;                   % 4    10    9    3
                    shape( 2)=   a*s1*s2*s3*t2*t3*t4;                   %
                    shape( 3)=   a*s1*s2*s3*t1*t2*t3;                   %
                    shape( 4)=   a*s2*s3*s4*t1*t2*t3;                   % 11   16   15    8
                    shape( 5)=-3.0*a*s1*s3*s4*t2*t3*t4;                 %
                    shape( 6)=-3.0*a*s1*s2*s4*t2*t3*t4;                 %
                    shape( 7)=-3.0*a*s1*s2*s3*t1*t3*t4;                 % 12   13   14    7
                    shape( 8)=-3.0*a*s1*s2*s3*t1*t2*t4;                 %
                    shape( 9)=-3.0*a*s1*s2*s4*t1*t2*t3;                 %
                    shape(10)=-3.0*a*s1*s3*s4*t1*t2*t3;                 % 1     5    6    2
                    shape(11)=-3.0*a*s2*s3*s4*t1*t2*t4;
                    shape(12)=-3.0*a*s2*s3*s4*t1*t3*t4;
                    shape(13)= 9.0*a*s1*s3*s4*t1*t3*t4;
                    shape(14)= 9.0*a*s1*s2*s4*t1*t3*t4;
                    shape(15)= 9.0*a*s1*s2*s4*t1*t2*t4;
                    shape(16)= 9.0*a*s1*s3*s4*t1*t2*t4;
                    deriv(1, 1)=  a *t2*t3*t4*(-s2*s3-s2*s4+s3*s4);
                    deriv(1, 2)=  a *t2*t3*t4*(-s1*s2+s1*s3+s2*s3);
                    deriv(1, 3)=  a *t1*t2*t3*(-s1*s2+s1*s3+s2*s3);
                    deriv(1, 4)=  a *t1*t2*t3*(-s2*s3-s2*s4+s3*s4);
                    deriv(1, 5)=-3.0*a*t2*t3*t4*(-s1*s3-s1*s4+s3*s4);
                    deriv(1, 6)=-3.0*a*t2*t3*t4*(-s1*s2+s1*s4+s2*s4);
                    deriv(1, 7)=-3.0*a*t1*t3*t4*(-s1*s2+s1*s3+s2*s3);
                    deriv(1, 8)=-3.0*a*t1*t2*t4*(-s1*s2+s1*s3+s2*s3);
                    deriv(1, 9)=-3.0*a*t1*t2*t3*(-s1*s2+s1*s4+s2*s4);
                    deriv(1,10)=-3.0*a*t1*t2*t3*(-s1*s3-s1*s4+s3*s4);
                    deriv(1,11)=-3.0*a*t1*t2*t4*(-s2*s3-s2*s4+s3*s4);
                    deriv(1,12)=-3.0*a*t1*t3*t4*(-s2*s3-s2*s4+s3*s4);
                    deriv(1,13)= 9.0*a*t1*t3*t4*(-s1*s3-s1*s4+s3*s4);
                    deriv(1,14)= 9.0*a*t1*t3*t4*(-s1*s2+s1*s4+s2*s4);
                    deriv(1,15)= 9.0*a*t1*t2*t4*(-s1*s2+s1*s4+s2*s4);
                    deriv(1,16)= 9.0*a*t1*t2*t4*(-s1*s3-s1*s4+s3*s4);
                    deriv(2, 1)=  a   *s2*s3*s4*(-t2*t3-t2*t4+t3*t4);
                    deriv(2, 2)=  a   *s1*s2*s3*(-t2*t3-t2*t4+t3*t4);
                    deriv(2, 3)=  a   *s1*s2*s3*(-t1*t2+t1*t3+t2*t3);
                    deriv(2, 4)=  a   *s2*s3*s4*(-t1*t2+t1*t3+t2*t3);
                    deriv(2, 5)= -3.0*a *s1*s3*s4*(-t2*t3-t2*t4+t3*t4);
                    deriv(2, 6)= -3.0*a *s1*s2*s4*(-t2*t3-t2*t4+t3*t4);
                    deriv(2, 7)= -3.0*a *s1*s2*s3*(-t1*t3-t1*t4+t3*t4);
                    deriv(2, 8)= -3.0*a *s1*s2*s3*(-t1*t2+t1*t4+t2*t4);
                    deriv(2, 9)= -3.0*a *s1*s2*s4*(-t1*t2+t1*t3+t2*t3);
                    deriv(2,10)= -3.0*a *s1*s3*s4*(-t1*t2+t1*t3+t2*t3);
                    deriv(2,11)= -3.0*a *s2*s3*s4*(-t1*t2+t1*t4+t2*t4);
                    deriv(2,12)= -3.0*a *s2*s3*s4*(-t1*t3-t1*t4+t3*t4);
                    deriv(2,13)=  9.0*a *s1*s3*s4*(-t1*t3-t1*t4+t3*t4);
                    deriv(2,14)=  9.0*a *s1*s2*s4*(-t1*t3-t1*t4+t3*t4);
                    deriv(2,15)=  9.0*a *s1*s2*s4*(-t1*t2+t1*t4+t2*t4);
                    deriv(2,16)=  9.0*a *s1*s3*s4*(-t1*t2+t1*t4+t2*t4);
                    if (neres==1)
                        heslo(1, 1)=    a *t2*t3*t4*(2.0*s2-2.0*s3-2.0*s4);
                        heslo(1, 2)=    a *t2*t3*t4*(-2.0*s1-2.0*s2+2.0*s3);
                        heslo(1, 3)=    a *t1*t2*t3*(-2.0*s1-2.0*s2+2.0*s3);
                        heslo(1, 4)=    a *t1*t2*t3*(2.0*s2-2.0*s3-2.0*s4);
                        heslo(1, 5)= -3.0*a *t2*t3*t4*(2.0*s1-2.0*s3-2.0*s4);
                        heslo(1, 6)= -3.0*a *t2*t3*t4*(-2.0*s1-2.0*s2+2.0*s4);
                        heslo(1, 7)= -3.0*a *t1*t3*t4*(-2.0*s1-2.0*s2+2.0*s3);
                        heslo(1, 8)= -3.0*a *t1*t2*t4*(-2.0*s1-2.0*s2+2.0*s3);
                        heslo(1, 9)= -3.0*a *t1*t2*t3*(-2.0*s1-2.0*s2+2.0*s4);
                        heslo(1,10)= -3.0*a *t1*t2*t3*(2.0*s1-2.0*s3-2.0*s4);
                        heslo(1,11)= -3.0*a *t1*t2*t4*(2.0*s2-2.0*s3-2.0*s4);
                        heslo(1,12)= -3.0*a *t1*t3*t4*(2.0*s2-2.0*s3-2.0*s4);
                        heslo(1,13)=  9.0*a *t1*t3*t4*(2.0*s1-2.0*s3-2.0*s4);
                        heslo(1,14)=  9.0*a *t1*t3*t4*(-2.0*s1-2.0*s2+2.0*s4);
                        heslo(1,15)=  9.0*a *t1*t2*t4*(-2.0*s1-2.0*s2+2.0*s4);
                        heslo(1,16)=  9.0*a *t1*t2*t4*(2.0*s1-2.0*s3-2.0*s4);
                        heslo(2, 1)=    a *s2*s3*s4*(2.0*t2-2.0*t3-2.0*t4);
                        heslo(2, 2)=    a *s1*s2*s3*(2.0*t2-2.0*t3-2.0*t4);
                        heslo(2, 3)=    a *s1*s2*s3*(-2.0*t1-2.0*t2+2.0*t3);
                        heslo(2, 4)=    a *s2*s3*s4*(-2.0*t1-2.0*t2+2.0*t3);
                        heslo(2, 5)= -3.0*a *s1*s3*s4*(2.0*t2-2.0*t3-2.0*t4);
                        heslo(2, 6)= -3.0*a *s1*s2*s4*(2.0*t2-2.0*t3-2.0*t4);
                        heslo(2, 7)= -3.0*a *s1*s2*s3*(2.0*t1-2.0*t3-2.0*t4);
                        heslo(2, 8)= -3.0*a *s1*s2*s3*(-2.0*t1-2.0*t2+2.0*t4);
                        heslo(2, 9)= -3.0*a *s1*s2*s4*(-2.0*t1-2.0*t2+2.0*t3);
                        heslo(2,10)= -3.0*a *s1*s3*s4*(-2.0*t1-2.0*t2+2.0*t3);
                        heslo(2,11)= -3.0*a *s2*s3*s4*(-2.0*t1-2.0*t2+2.0*t4);
                        heslo(2,12)= -3.0*a *s2*s3*s4*(2.0*t1-2.0*t3-2.0*t4);
                        heslo(2,13)=  9.0*a *s1*s3*s4*(2.0*t1-2.0*t3-2.0*t4);
                        heslo(2,14)=  9.0*a *s1*s2*s4*(2.0*t1-2.0*t3-2.0*t4);
                        heslo(2,15)=  9.0*a *s1*s2*s4*(-2.0*t1-2.0*t2+2.0*t4);
                        heslo(2,16)=  9.0*a *s1*s3*s4*(-2.0*t1-2.0*t2+2.0*t4);
                        heslo(3, 1)=    a*(-s2*s3-s2*s4+s3*s4)*(-t2*t3-t2*t4+t3*t4);
                        heslo(3, 2)=    a*(-s1*s2+s1*s3+s2*s3)*(-t2*t3-t2*t4+t3*t4);
                        heslo(3, 3)=    a*(-s1*s2+s1*s3+s2*s3)*(-t1*t2+t1*t3+t2*t3);
                        heslo(3, 4)=    a*(-s2*s3-s2*s4+s3*s4)*(-t1*t2+t1*t3+t2*t3);
                        heslo(3, 5)= -3.0*a*(-s1*s3-s1*s4+s3*s4)*(-t2*t3-t2*t4+t3*t4);
                        heslo(3, 6)= -3.0*a*(-s1*s2+s1*s4+s2*s4)*(-t2*t3-t2*t4+t3*t4);
                        heslo(3, 7)= -3.0*a*(-s1*s2+s1*s3+s2*s3)*(-t1*t3-t1*t4+t3*t4);
                        heslo(3, 8)= -3.0*a*(-s1*s2+s1*s3+s2*s3)*(-t1*t2+t1*t4+t2*t4);
                        heslo(3, 9)= -3.0*a*(-s1*s2+s1*s4+s2*s4)*(-t1*t2+t1*t3+t2*t3);
                        heslo(3,10)= -3.0*a*(-s1*s3-s1*s4+s3*s4)*(-t1*t2+t1*t3+t2*t3);
                        heslo(3,11)= -3.0*a*(-s2*s3-s2*s4+s3*s4)*(-t1*t2+t1*t4+t2*t4);
                        heslo(3,12)= -3.0*a*(-s2*s3-s2*s4+s3*s4)*(-t1*t3-t1*t4+t3*t4);
                        heslo(3,13)=  9.0*a*(-s1*s3-s1*s4+s3*s4)*(-t1*t3-t1*t4+t3*t4);
                        heslo(3,14)=  9.0*a*(-s1*s2+s1*s4+s2*s4)*(-t1*t3-t1*t4+t3*t4);
                        heslo(3,15)=  9.0*a*(-s1*s2+s1*s4+s2*s4)*(-t1*t2+t1*t4+t2*t4);
                        heslo(3,16)=  9.0*a*(-s1*s3-s1*s4+s3*s4)*(-t1*t2+t1*t4+t2*t4);
                    end
                end
        end
    case '3D'
        ndime = 3;
        s=posgp(1,igaus);t=posgp(2,igaus);u=posgp(3,igaus);
        deriv = zeros(ndime,nnode);
        shape = zeros(1,nnode);
        heslo = zeros(3,nnode);
        switch  etype
            case 'HEXAHEDRA'
                if(nnode==8)
                    lcord(1,1)= -1; lcord(1,2)= -1; lcord(1,3)= -1;
                    lcord(2,1)=  1; lcord(2,2)= -1; lcord(2,3)= -1;
                    lcord(3,1)=  1; lcord(3,2)=  1; lcord(3,3)= -1;
                    lcord(4,1)= -1; lcord(4,2)=  1; lcord(4,3)= -1;
                    lcord(5,1)= -1; lcord(5,2)= -1; lcord(5,3)=  1;
                    lcord(6,1)=  1; lcord(6,2)= -1; lcord(6,3)=  1;
                    lcord(7,1)=  1; lcord(7,2)=  1; lcord(7,3)=  1;
                    lcord(8,1)= -1; lcord(8,2)=  1; lcord(8,3)=  1;
                    for inode=1:nnode
                        shape(inode)=(1+lcord(inode,1)*s)*(1+lcord(inode,2)*t)*(1+lcord(inode,3)*u)/8;
                        deriv(1,inode)=lcord(inode,1)*(1+lcord(inode,2)*t)*(1+lcord(inode,3)*u)/8;
                        deriv(2,inode)=lcord(inode,2)*(1+lcord(inode,1)*s)*(1+lcord(inode,3)*u)/8;
                        deriv(3,inode)=lcord(inode,3)*(1+lcord(inode,1)*s)*(1+lcord(inode,2)*t)/8;
                    end
                end
        end
        
end

end

