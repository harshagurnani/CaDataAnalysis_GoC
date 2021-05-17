function newParams = basicplot_params( oldParams )
%% Params for BasicPlot - Traces and Patches
% Add defaults where params are missing

% ----------------
BasicPlot.two_col          = true;
BasicPlot.scale            = 30;                                       % 0 if it should be determined automatically
BasicPlot.showPatch        = true;
BasicPlot.graylim          = [0 500];                                  % 0 if it should be determined automatically - for gray patch mean image
BasicPlot.sortDepth        = 'descend';                                % Descend means deepest neurons are lowest traces. 'ascend' is reverse. If empty. then it doesn't sort.
BasicPlot.addPuff          = true;                                     
BasicPlot.addSpeed         = true;

BasicPlot.LineWidth        = 1.5;
BasicPlot.ColorMap         = 'jet';                                    % Can be a matrix
BasicPlot.EventColor.Loco  = 'k';
BasicPlot.EventColor.Whisk = 'g';
BasicPlot.EventColor.Puff  = 'm';

BasicPlot.dff_scale        = 30;                          % factor to multiply dff 
BasicPlot.graylim          = [50 800];                    % mean patch image - range
BasicPlot.cropTBins        = 10;                          % how many timepoints to drop from the beginning?
BasicPlot.TmOffset         = 2;                           % second% shift traces along x-axis (space between patch and dff trace)


BasicPlot.patchSize        = [10 50];
% merge defaults into struct
allfields = fieldnames( BasicPlot);
for ff = 1:numel(allfields )
   if ~isfield( oldParams, allfields{ff} )
       oldParams.(allfields{ff}) = BasicPlot.(allfields{ff}); 
   end
end

newParams = oldParams;
end