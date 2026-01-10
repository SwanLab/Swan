function [PhiMaster_nonl,PhiSlave_nonl,qMASTER_nonl,...
    qSLAVE_nonl,...
    DATA_evaluateTAU_and_DER, nREDcoor] = Bsplines_EP_1stmodeIND(Phi_non,Kll,DATA_interp,SNAPdisp_plast_L)

if nargin == 0
    load('tmp1.mat')
    close all
    DATA_interp.order_Bslines = 4; 
end
DATA_interp = DefaultField(DATA_interp, 'order_Bslines', 4);

order_Bslines = DATA_interp.order_Bslines;

% % This is the master mode (the first mode)
% PhiMaster_nonl =    Phi_non(:,1) ;
% % Remaining modes (slave modes)
% PhiSlave_nonl =  Phi_non(:,2:end) ; % Slave modes, plastic
% % Slave coefficients 
% qSLAVE_nonl = PhiSlave_nonl'*(Kll*SNAPdisp_plast_L) ;
% % Master coefficients
% qMASTER_nonl = PhiMaster_nonl'*(Kll*SNAPdisp_plast_L) ;


% IN THIS "NEW" APPROACH, WE FIRST CONDUCT A REGRESSION PROBLEM
% IN WHICH WE APPROXIMATE qMASTER_nonl and qSLAVE_nonl 
% as a function of the input parameter
disp(['Pre-fitting (only valid for case in which  qSLAVE_nonl is monotonic ...)'])

q_nonl = Phi_non'*(Kll*SNAPdisp_plast_L) ;

q_REFERENCE = linspace(0,1,size(q_nonl,2)) ; 
  figure(200  )
  hold on 
  title(['Original qNON and Bspline approximations '])
  xlabel('q REFERENCE')
  ylabel('q')
  
  q_nonlBSPL = q_nonl ; 
  
  sp = cell(size(q_nonl,1),1) ; 
knots_number = DATA_interp.NSAMPLES ;
colors = lines(size(q_nonl,1));


for islave = 1:size(q_nonl,1)
    sp{islave} = spap2(knots_number, order_Bslines, q_REFERENCE, q_nonl(islave,:)); % Function
    sp1 =  fnder(sp{islave}) ;  % First derivative
    sp2 =  fnder(sp1) ;  % Second derivative
    
    
    [OUTPUT ] = evaluate_spline_with_extrapolation( sp{islave} , sp1 ,  sp2, q_REFERENCE) ;
    q_nonlBSPL(islave,:) = OUTPUT.VALUE ; 
   fff= max(abs(q_nonl(islave,:))); 
   
   q_nonlBSPL(islave,:) = OUTPUT.VALUE ; 
 
    plot(q_REFERENCE, q_nonl(islave,:)/fff,'DisplayName',['q',num2str(islave),' MAXq=',num2str(fff)],'Color',colors(islave,:)) ;
    
    plot(q_REFERENCE,OUTPUT.VALUE/fff,'--','DisplayName',['Bspline ',num2str(islave)],'Color',colors(islave,:)) ;
   
    
end

legend show


% NEXT WE SELECT THE "MASTER" mode


% This is the master mode (the first mode)
PhiMaster_nonl =    Phi_non(:,1) ;
% Remaining modes (slave modes)
PhiSlave_nonl =  Phi_non(:,2:end) ; % Slave modes, plastic
% Slave coefficients 
qSLAVE_nonl = q_nonlBSPL(2:end,:) ;
% Master coefficients
qMASTER_nonl =  q_nonlBSPL(1,:) ; 


% 
% 
% figure(35)
% hold on
% xlabel('Snapshots')
% ylabel('q')
% qALL = [qMASTER_nonl;qSLAVE_nonl] ;
% for iii  =1:size(qALL,1)
%     plot( qALL(iii,:),'DisplayName',['q',num2str(iii)])
% end
% legend show

% B-splines, it is necessary to sort the entries of qMASTER_nonl is
% ascending order
%         [qMASTER_nonl,indORDER] = sort(qMASTER_nonl) ;
%
%         figure(126)
%         hold on
%         xlabel('Snapshot ordered')
%         ylabel('Master function (qMASTER nonl)')
%         plot(qMASTER_nonl)
%
%         qSLAVE_nonl = qSLAVE_nonl(:,indORDER)  ;
%
% Now we have to establish the mapping qSLAVE_nonl = f(qMASTER_nonl)
% We do this in an indirect way, by first making
%% qSLAVE_nonl = UU*diag(SS)*VV^T
% and then constructing the mapping VV(:,i) = f(qMASTER_nonl)
%
%

% SAMPLING RATE 
qMASTER_orig = q_nonl(1,:) ; 
[qSORT_master,IndSort] = sort(qMASTER_orig) ; 
Inc_qSORT = diff(qSORT_master) ; 
FACTOR = 3; 
SamplingRate =  max(Inc_qSORT)/FACTOR ; 
NstepsLOC = ceil((qSORT_master(end)-qSORT_master(1))/SamplingRate) ; 
SAmpling_qMASTER = linspace(qSORT_master(1),qSORT_master(end),NstepsLOC) ; 
[AAA,BBB] = knnsearch(qSORT_master',SAmpling_qMASTER') ; 

IndSelect = unique(AAA) ; 
 IndSort_new = IndSort(IndSelect) ; 
 
 
qSLAVE_nonl = q_nonlBSPL(2:end,IndSort_new) ;
% Master coefficients
qMASTER_nonl =  q_nonlBSPL(1,IndSort_new) ; 




[UU,SS,VV] = SVDT(qSLAVE_nonl) ;

if length(qSLAVE_nonl) >= DATA_interp.NSAMPLES
    DATA_interp.NSAMPLES = ceil(0.9*length(qSLAVE_nonl)) ; 
end

% figure(345)
% hold on
% title('Right singular vectors qSLAVE nonl as a function of qMASTER nonl')
% xlabel('qMASTER nonl')
% ylabel('qSLAVE nonl')
%
%
% for iii = 1:size(VV,2)
%     plot(qMASTER_nonl,VV(:,iii),'DisplayName',['\lambda',num2str(iii),'=',num2str(SS(iii)/SS(1))]) ;
% end


[DATA_evaluateTAU_and_DER, nREDcoor] = BsplinesLeastSquares_PLAST(DATA_interp, qMASTER_nonl, VV', UU, SS) ;
