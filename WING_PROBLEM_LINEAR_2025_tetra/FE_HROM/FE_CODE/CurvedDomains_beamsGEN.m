function   DATA3D  = CurvedDomains_beamsGEN(DATA3D,angDOM,DATA)

if nargin == 0
    load('tmp2.mat')
end
DATA = DefaultField(DATA,'ISTWIST_ANGLE',0) ;   
if DATA.ISTWIST_ANGLE == 0
    DATA3D  = CurvedDomains_beams(DATA3D,angDOM) ; 
else
    DATA3D  = CurvedDomains_beamsTWIST(DATA3D,angDOM) ; 
end
 