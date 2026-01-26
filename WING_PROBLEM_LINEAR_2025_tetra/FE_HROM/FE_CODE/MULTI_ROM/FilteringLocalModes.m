function [TEXTB,MODES_INCLUDE,IMODES_INCLUDE,ratioSV_ABS_1,ratioSV_ABS_2] = FilteringLocalModes(ratioSV_f1,TOL_SINGULAR_VALUES_Hqr,ratioSV_f2,SSVAL_f1,SSVAL_f2,sNEWMODES,TOL_REL_MAX_SV,...
    NEW_MODES,imode,IMODES_INCLUDE,TEXTB,MODES_INCLUDE,DATAIN,ratioSV_ABS_1,ratioSV_ABS_2)

%if DATAIN.FILTERING_INF_MODES_MINIMUM_CORRELATION_ON_JUST_ONE_FACE  == 0
ratioSV_ABS_1 = 0 ; 
ratioSV_ABS_2 = 0 ; 
if (ratioSV_f1 >= TOL_SINGULAR_VALUES_Hqr && ratioSV_f2 >= TOL_SINGULAR_VALUES_Hqr) ...
        && length(SSVAL_f1) == sNEWMODES && length(SSVAL_f2) == sNEWMODES
    if isempty(TOL_REL_MAX_SV)
        MODES_INCLUDE = NEW_MODES ;
        IMODES_INCLUDE(end+1) = imode ;
    else
        ratioSV_ABS_1 = SSVAL_f1(end)/SSVAL_f1(1) ;
        ratioSV_ABS_2 = SSVAL_f2(end)/SSVAL_f2(1) ;
        if ratioSV_ABS_1 >= TOL_REL_MAX_SV && ratioSV_ABS_2 >= TOL_REL_MAX_SV
            MODES_INCLUDE = NEW_MODES ;
            IMODES_INCLUDE(end+1) = imode ;
        else
            TEXTB{end+1} = ['Mode ',num2str(imode),' excluded because S(end)/S(1) exceeds the tolerance'];
            TEXTB{end+1} = ['[S(end)/S(1)]_f1 =',num2str(ratioSV_ABS_1),'; ','[S(end)/S(1)]_f2 =',num2str(ratioSV_ABS_2)];
        end
    end
else
    
    TEXTB{end+1} = ['Mode ',num2str(imode),' excluded because S(end-1)/S(end) exceeds the tolerance'];
    TEXTB{end+1} = ['[S(end)/S(end-1)]_f1 =',num2str(ratioSV_f1),'; ','[S(end)/S(end-1)]_f2 =',num2str(ratioSV_f2)]
end
%
% else
%     %%% Change && for ||
%
%       if (ratioSV_f1 >= TOL_SINGULAR_VALUES_Hqr || ratioSV_f2 >= TOL_SINGULAR_VALUES_Hqr) ...
%             && length(SSVAL_f1) == sNEWMODES || length(SSVAL_f2) == sNEWMODES
%         if isempty(TOL_REL_MAX_SV)
%             MODES_INCLUDE = NEW_MODES ;
%             IMODES_INCLUDE(end+1) = imode ;
%         else
%             ratioSV_ABS_1 = SSVAL_f1(end)/SSVAL_f1(1) ;
%             ratioSV_ABS_2 = SSVAL_f2(end)/SSVAL_f2(1) ;
%             if ratioSV_ABS_1 >= TOL_REL_MAX_SV || ratioSV_ABS_2 >= TOL_REL_MAX_SV
%                 MODES_INCLUDE = NEW_MODES ;
%                 IMODES_INCLUDE(end+1) = imode ;
%             else
%                 TEXTB{end+1} = ['Mode ',num2str(imode),' excluded because S(end)/S(1) exceeds the tolerance'];
%                 TEXTB{end+1} = ['[S(end)/S(1)]_f1 =',num2str(ratioSV_ABS_1),'; ','[S(end)/S(1)]_f2 =',num2str(ratioSV_ABS_2)];
%             end
%         end
%     else
%
%         TEXTB{end+1} = ['Mode ',num2str(imode),' excluded because S(end-1)/S(end) exceeds the tolerance'];
%         TEXTB{end+1} = ['[S(end)/S(end-1)]_f1 =',num2str(ratioSV_f1),'; ','[S(end)/S(end-1)]_f2 =',num2str(ratioSV_f2)]
%     end
%
%
% end