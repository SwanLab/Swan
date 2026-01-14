function [DATA_evaluateTAU_and_DER,DATA_interp] =  BsplinesSUBSAMPL_damage(qMASTER_nonl,DATA_interp,qSLAVE_nonl)
%  BsplinesSUBSAMPL_damage is a modification of   BsplinesSUBSAMPL
% to deal with data generated via a damage model
%  JAHO, 23-Oct-2025, Thursday,  Balmes 185, Barcelona.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/17_DAMAGE.mlx
if nargin == 0
    load('tmp1.mat')
end

% Sampling uniform points from the image of  qMASTER_nonl
qMASTER_range = linspace(min(qMASTER_nonl),max(qMASTER_nonl),DATA_interp.NSAMPLES) ;
[IndSelect,~] = knnsearch(qMASTER_nonl',qMASTER_range') ;
IndSelectUNI = unique(IndSelect) ;

qMASTER_nonl_SAMP = qMASTER_nonl(IndSelectUNI) ;
qSLAVE_nonl_SAMP = qSLAVE_nonl(:,IndSelectUNI) ;


figure(35)
subplot(2,1,1)

hold on
xlabel('Snapshots')
ylabel('Modal amplitudes')
title(['Assessing the effect of subsampling (nsamples = ',num2str(DATA_interp.NSAMPLES),')'])


subplot(2,1,2)
hold on
xlabel('Snapshots')
ylabel('Der. Modal amplitudes')

qALL = [qMASTER_nonl;qSLAVE_nonl] ;
qSUB = [qMASTER_nonl_SAMP;qSLAVE_nonl_SAMP] ;
colors = lines(size(qSUB,1));

DATA_interp= DefaultField(DATA_interp,'ratio_NSAMPLES_knots',1);

knots_max=  length(IndSelectUNI)-DATA_interp.order_Bslines+1;

knots_max = knots_max*DATA_interp.ratio_NSAMPLES_knots ;

%n_knots = ceil(DATA_interp.ratio_NSAMPLES_knots*length(IndSelectUNI)) ;

for iii  =1:size(qALL,1)
    subplot(2,1,1)
    hold on
    plot( qALL(iii,:),'DisplayName',['qORIG',num2str(iii)],'Color',colors(iii,:)) ;
    
    sp = spap2(knots_max, DATA_interp.order_Bslines, IndSelectUNI, qSUB(iii,:)); % Function
    sp1 =  fnder(sp ) ;  % First derivative
    sp2 =  fnder(sp1) ;  % Second derivative
    
    
    [OUTPUT ] = evaluate_spline_with_extrapolation( sp, sp1 ,  sp2, 1:size(qALL,2)) ;
    %     q_nonlBSPL(islave,:) = OUTPUT.VALUE ;
    %    fff= max(abs(q_nonl(islave,:)));
    %
    %    q_nonlBSPL(islave,:) = OUTPUT.VALUE ;
    %
    
    
    plot(1:size(qALL,2), OUTPUT.VALUE,'--','DisplayName',['qSUB splin',num2str(iii)],'Color',colors(iii,:)) ;
    
    subplot(2,1,2)
    hold on
    plot(1:size(qALL,2), OUTPUT.DER1,'--','DisplayName',['qSUB splin',num2str(iii)],'Color',colors(iii,:)) ;
    
end
legend show

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

%


[UU,SS,VV] = SVDT(qSLAVE_nonl_SAMP) ;

%   DATA_interp= DefaultField(DATA_interp,'ratio_NSAMPLES_knots',0.9);

DATA_interp.NSAMPLES =knots_max;
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
[DATA_evaluateTAU_and_DER] = BsplinesLeastSquares_DAMAGE(DATA_interp, qMASTER_nonl_SAMP, VV', UU, SS) ;