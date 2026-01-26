function uBARloc = PrescribedRIGIBbodyDisp_3D(BasisUrb,uBAR_0) 


 
    uBARloc = zeros(length(FACE),1) ;
    uBARloc(1:3:end) = uBAR_0(1) ;
    uBARloc(2:3:end) = uBAR_0(2) ;
    uBARloc(3:3:end) = uBAR_0(3) ;
    uBARloc= uBARloc + uBAR_0(4)*BasisUrb(FACE,4) ;
    uBARloc= uBARloc + uBAR_0(5)*BasisUrb(FACE,5) ;
    uBARloc= uBARloc + uBAR_0(6)*BasisUrb(FACE,6) ;
 