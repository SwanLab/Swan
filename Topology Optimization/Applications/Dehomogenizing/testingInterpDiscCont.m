function testingInterpDiscCont    
close all
connec = [1 2 4;
          1 4 3];
coord = [0 0;
         1 0;
         0 1;
         1 1];
sC.connec = connec;
sC.coord  = coord;
mC = Mesh(sC);


beta = pi/180*[190;175;60;0];
alpha = beta/2;
tC = [cos(alpha) sin(alpha)];

plotOrientation(mC,tC)


sD.connec = reshape(1:mC.nnodeElem*mC.nelem,mC.nnodeElem,mC.nelem)';
d = connec';
sD.coord(:,1) = coord(d(:),1);
sD.coord(:,2) = coord(d(:),2);
mD = Mesh(sD);

tD(:,1) = tC(d(:),1);
tD(:,2) = tC(d(:),2);

plotOrientation(mD,tD)

s.meshCont   = mC;
s.meshDisc   = mD;
s.fieldCont  = tC;
s.fieldDisc  = tD;

c = ContDiscontinousInterpolator(s);
c.compute();

end

function plotOrientation(m,t)
figure()
m.plot()
quiver(m.coord(:,1),m.coord(:,2),t(:,1),t(:,2))
end