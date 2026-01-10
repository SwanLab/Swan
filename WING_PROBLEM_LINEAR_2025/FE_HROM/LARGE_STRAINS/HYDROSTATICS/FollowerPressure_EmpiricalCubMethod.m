function [ECMdata,Basis_tPRESS,S_tPRESS] = FollowerPressure_EmpiricalCubMethod(NAMEsnap_base,DATAoffline,BasisU,CASES,DOFr,DOFl,DATAparamSTUDY,OPERFE,MATPRO,DATA)
% Counterpart of  Stress__EmpiricalCubMethod for follower loads
% JAH0, 5-July-2021 ---- UPCT, Cartagena
% -----------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end


% Follower pressure loads as a function of velocities and displacements
% ------------------------------------------------------------------------------
[Basis_tPRESS,S_tPRESS] = tPRESS_fromBasisFromDisplacements(BasisU,CASES,NAMEsnap_base,DOFr,DOFl,DATAparamSTUDY,OPERFE,MATPRO,DATAoffline,DATA) ;
% -------------------
% HYPERREDUCTION
% -------------------
%
% else

%
NstRED_l = ConvertBlockDiag_general(OPERFE.HYDRO.NbST,DATA.MESH.ndim)*OPERFE.HYDRO.Lbool(:,DOFl)*BasisU ;

%
% end
%if DATAoffline.USE_ELEMENT_BASED_APPROACH == 0
% ******************
% Discrete ECM
% ******************
%   DATA_GENGAUSS.ACTIVE =  0 ;
%   if DATA_GENGAUSS.ACTIVE == 0
if ~isempty(Basis_tPRESS)
[setPoints,wRED,ERROR_GLO,DATAOUT] = DiscreteECM_follower(NstRED_l,Basis_tPRESS,DATA,OPERFE.HYDRO.wST,DATAoffline) ;
ECMdata.setPoints = setPoints ;
ECMdata.wRED = wRED ;
proporP = length(setPoints)/DATA.MESH.HYDRO.ngausT*100;
disp(['Number of ECM points = ',num2str(length(setPoints)),' (',num2str(proporP),' % total)'])

setElements = large2smallREP(setPoints,DATA.MESH.HYDRO.ngausT/DATA.MESH.HYDRO.nelemB) ;

setElementsSHOW = setElements + DATA.MESH.nelem ;  
disp('****************************+')
disp(['List of boundary elements  = ',num2str(length(setElementsSHOW)),' elements'])
disp(num2str(setElementsSHOW'))
%clipboard('copy',num2str(setElements'));
ECMdata.setElements = setElements ;

else
   ECMdata = [] ; 
    
end

 