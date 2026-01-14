function Se = ChangeCoordBnd(Xe,ndimB)
% Given the   matrix (Xe) of global coordinates of the position of nodes of boundary element "e", this function
% expresses these positions in a coordinate system intrinsic to the
% element (Se).
if nargin == 0
    load('tmp2.mat')
end
nnodeEb = size(Xe,2) ;
t1 = Xe(:,2)-Xe(:,1) ;    t1 = t1/norm(t1) ;
R = zeros(ndimB+1,ndimB+1) ;
if ndimB==1
    R(:,1) = t1 ; R(:,2)=[-t1(2);t1(1)];
elseif ndimB == 2
    R(:,1) = t1 ;
    t2 = Xe(:,3)-Xe(:,1) ;    t2 = t2/norm(t2) ;
    t1t2 = cross(t1,t2);
    R(:,3) =  t1t2/norm(t1t2) ;
    R(:,2) = cross(R(:,3),R(:,1)) ;
end
 Se = zeros(ndimB,nnodeEb) ; 
    for i=1:nnodeEb
        xeR = R'*(Xe(:,i)-Xe(:,1)) ;
        Se(:,i) = xeR(1:ndimB) ; 
    end   