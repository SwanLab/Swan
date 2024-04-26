function intersect=TriangleTriangleIntersection(P1,P2,P3,O1,O2,O3,ignore_corners)
for i=1:6
    switch(i)
        case 1
            A=P1; B=P2; C=P3; E1=O1; E2=O2;
        case 2
            A=P1; B=P2; C=P3; E1=O2; E2=O3;
        case 3
            A=P1; B=P2; C=P3; E1=O3; E2=O1;
        case 4
            A=O1; B=O2; C=O3; E1=P1; E2=P2;
        case 5
            A=O1; B=O2; C=O3; E1=P2; E2=P3;
        case 6
            A=O1; B=O2; C=O3; E1=P3; E2=P1;
    end
    intersect=LineTriangleIntersection(A,B,C,E1,E2,ignore_corners);
    if(intersect), return; end
end
