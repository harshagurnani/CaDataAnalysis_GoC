function [ all_traces, time_trace ] = avg_grouped_roi( group_struct, raw, times, params, varargin)
%% To group traces in same soma or dendrite, and reorder dendrites next to somas they belong to. Floating parts are added at the end.
%
% ARGUMENTS - 
%   - group_struct      - Structure with 3 cell arrays - soma, dend, and
%                           other_groups. Each cell in soma is numeric array of
%                           ROI numbers belonging to the soma. Each cell of
%                           dend is a cell array of dendrite ROIs belonging to
%                           the soma with same array index. Each cell in the
%                           soma-specific cell array is numeric array of points
%                           on a single branch to be averaged. Indexing:
%                           group.soma{1} - array of points to be averaged for 1st soma
%                           group.dend{1}{3} - array of points to be averaged for 3rd dendritic segment/branch of 1st soma
%                           group.other_groups{2} - 2nd floating segment
%   - raw               - Raw data (3-dim array : Timepoints x Trials x ROI)
%   - times             - Time data (3-dim array : Timepoints x Trials x ROI)
%   - params            - Params structure (modified or same as wfor loading data
%   
%   OPTIONAL
%   - donorm (LOGICAL)      - Normalise group-averaged traces
%   - doplot (LOGICAL)      - Plot averaged (and normalised) traces
%   - zscore_plot (LOGICAL) - Z-score agroup-averaged traces and plot heatmap-like plot
%   - concat (LOGICAL)      - Concatenate all_traces into a single matrix
%
% RETURNS:
%   - all_traces (Cell array or 3D array) - Averaged traces for each cell (or other
%                               group). Somatic and dendritic traces of the
%                               same marked cell are in the same cell array
%                               if concat=FALSE. Otherwise a single 3D
%                               array with diff groups as diff ROIs is
%                               returned.
%   - time_trace (Cell array or 3D array) - Time vector for each averaged group
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prevn=0;
nsoma =  size(group_struct.soma,2);
if isfield( group_struct, 'other_groups')
    ngroups =  size(group_struct.other_groups,2);
else
    ngroups = 0;
end
[ donorm, doplot, zscore_plot, concat ] = parse_optional_args( varargin, nargin, 4);

if doplot
    figure;
end
for celln=1:nsoma
    s_trace = mean(raw(:,:,group_struct.soma{celln}(:)),3);
    d_trace = [];
    ndend = size(group_struct.dend{celln},2);
    for dd = 1:ndend
        d_trace = cat( 3, d_trace, mean(raw(:,:,group_struct.dend{celln}{dd}),3));
    end
    all_traces{celln} = cat(3, s_trace, d_trace);

    time_trace{celln} = mean(times(:,:,group_struct.soma{celln}(:)),3);
    for dd = 1:ndend
        time_trace{celln} = cat( 3, time_trace{celln}, mean(times(:,:,group_struct.dend{celln}{dd}),3));
    end
    
    if donorm
        all_traces{celln} = normalize_and_smooth( all_traces{celln}, time_trace{celln}, params);
    end
    if doplot
        plot_traces( all_traces{celln}+prevn+1*(celln-1), time_trace{celln}, [ 1 3 2])
        hold on
    end
    prevn = prevn+1+ndend;
    
    
end

for celln=1:ngroups
    if iscell(group_struct.other_groups{celln})
        ndend = size(group_struct.other_groups{celln},2);
        d_trace = [];
        time_trace{celln+nsoma} = [];
        for dd = 1:ndend
            d_trace = cat( 3, d_trace, mean(raw(:,:,group_struct.other_groups{celln}{dd}),3));
            time_trace{celln+nsoma} = cat( 3, time_trace{celln+nsoma}, mean(times(:,:,group_struct.other_groups{celln}{dd}),3));
        end
        all_traces{celln+nsoma} =  d_trace;    
    else
        ndend = 0;
        all_traces{celln+nsoma} = mean(raw(:,:,group_struct.other_groups{celln}(:)),3);
        time_trace{celln+nsoma} = mean(times(:,:,group_struct.other_groups{celln}(:)),3);
    end
    if donorm
        all_traces{celln+nsoma} = normalize_and_smooth( all_traces{celln+nsoma}, time_trace{celln+nsoma}, params);
    end
    if doplot
        plot_traces( all_traces{celln+nsoma}+prevn+1*(celln+nsoma-1), time_trace{celln+nsoma}, [ 1 3 2])
        hold on
    end
    prevn = prevn+1+ndend;
end


 if concat || zscore_plot
     all_concat = [];
     for celln = 1:nsoma+ngroups
         all_concat = cat( 3, all_concat, all_traces{celln} );
     end
 end
    
 if concat
     all_traces = all_concat;
 end
 
 if zscore_plot
     figure;
     plot_scaled_image( all_concat, [1 3 2], 'zscore')
 end

end

%% Parse optional arguments
function [ donorm, doplot, zscore_plot, concat ] = parse_optional_args( varargs, nargs, copts)

donorm=[];
doplot = [];
zscore_plot = [];
concat = [];

if nargs > copts  

    charctr=0;
    numctr=0;
    start_read_var = copts+1;
    while start_read_var <= nargs
        vn = start_read_var-copts;
       if ischar( varargs{vn} )
           charctr=1;
           if strcmp(varargs{vn}, 'donorm' )
               donorm = varargs{vn+1};
               start_read_var = start_read_var + 2;
           elseif strcmp(varargs{vn}, 'doplot' )
               doplot = varargs{vn+1};
               start_read_var = start_read_var + 2;  
           elseif strcmp(varargs{vn}, 'zscore_plot' )
               zscore_plot = varargs{vn+1};
               start_read_var = start_read_var + 2;  
           elseif strcmp(varargs{vn}, 'concat' )
               concat = varargs{vn+1};
               start_read_var = start_read_var + 2;  
           else 
                error('Argument name %s is not an optional argument', varargs{vn})
           end
       else
           if charctr == 1
               error('Paired name-value cannot come before value arguments')
           end
           numctr = numctr+1;
           switch numctr
               case 1
                   donorm = varargs{vn};
                   start_read_var = start_read_var+1;
               case 2
                   doplot = varargs{vn};
                   start_read_var = start_read_var+1;
               case 3
                   zscore_plot = varargs{vn};
                   start_read_var = start_read_var+1;
               case 4
                   concat = varargs{vn};
                   start_read_var = start_read_var+1;
           end
       end
    end

end

if isempty(donorm)
    donorm = true;
end

if isempty(doplot)
    doplot=true;
end

if isempty(zscore_plot)
    zscore_plot=false;
end

if isempty(concat)
    concat=false;
end
end