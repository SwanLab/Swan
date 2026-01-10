function ECMdata = SAW_ECM_large1param_CONT(DATAoffline,SNAPfint,DATA,OPERFE,qLATENT,DATA_interp)
% JAHO, 13-Sept-2025, Saturday, Thriller, Caffe bar, Sarajevo, Bosnia.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/10_SAW_ECMlarge.mlx
if nargin == 0
    load('tmp3.mat')
    %   DATAoffline.IntegrateExactlyVolume_SAW_ECM = 1;
    close all
end


DATAoffline= DefaultField(DATAoffline,'UsePreliminaryVersion_SAW_ECM',0) ;


if  DATAoffline.UsePreliminaryVersion_SAW_ECM == 1
    % Not adaptive per se yet
    [ECMdata_cluster,setCandidates,qLATENT,SNAPfint] = SAW_ECM_large1param_CONT_prelim(DATAoffline,SNAPfint,DATA,OPERFE,qLATENT) ;
    
    
    disp(['Number of Candidate Points'])
    disp(length(setCandidates))
    disp(['setCandidates = '])
    disp(setCandidates')
    
    
    %%%%%%%%%%%%%%%%%5
    %Creating the chart qLATENT-wECM,zECM
    setPointsALL = cell(size(ECMdata_cluster)) ;
    
    for icluster = 1:length(ECMdata_cluster)
        setPointsALL{icluster} =  ECMdata_cluster{icluster}.setPoints  ;
    end
    setPointsALL  =unique(cell2mat(setPointsALL' )) ;
    
    ncluster = length(ECMdata_cluster) ;
    npointsALL = length(setPointsALL) ;
    wALL = zeros(npointsALL,ncluster) ;
    
    
    
    for icluster = 1:ncluster
        setPloc = ECMdata_cluster{icluster}.setPoints  ;
        [dummy1,III,JJJ] = intersect(setPloc,setPointsALL,'stable') ;
        wALL(JJJ,icluster) =  ECMdata_cluster{icluster}.wRED  ;
        
        
    end
    
    figure(4532)
    hold on
    title('WEIGHTS VERSUS LATENT VARIABLE (HYPERELASTIC PROBLEM)')
    xlabel('q')
    ylabel('w')
    for iddd  = 1:size(wALL,1)
        eeee = large2smallINCLUDEREP(setPointsALL(iddd),DATA.MESH.ngaus) ;
        plot(qLATENT,wALL(iddd,:),'DisplayName',['Point =',num2str(setPointsALL(iddd)),' Elem = ',num2str(eeee)])
    end
    legend show
    
    
    
    
    
    
    
    
    
else
    % Adaptive ECM (MAW-ECM)
    [wALL,setPointsALL,qLATENT,SNAPfint] = ManifoldAdaptWeightsECM(DATAoffline,SNAPfint,DATA,OPERFE,qLATENT) ;
end



% CONSTRUCT THE REGRESSION MAPPING FROM qLATENT to wALL
DATA_regress = MAW_ECM_regression(qLATENT,wALL,DATA_interp ) ;


ECMdata.wRED.DATA_regress = DATA_regress;
ECMdata.setPoints = setPointsALL ;
% ECMdata.wRED.Values = wALL ;
% ECMdata.wRED.q = qLATENT ;
setElements_cand = large2smallINCLUDEREP(setPointsALL,DATA.MESH.ngaus) ;
ECMdata.setElements = setElements_cand ;
ECMdata.wRED.IndexDOFl_q = 1;

disp('Cheching accuracy MAW-ECM')
CheckAccurayMAT_ECM(ECMdata,SNAPfint,OPERFE,qLATENT) ;

