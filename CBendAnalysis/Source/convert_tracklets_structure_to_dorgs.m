% This code is part of the MATLAB toolbox Gait-CAD.
% Copyright (C) 2012 [Johannes Stegmaier, Ralf Mikut]
%
%
% Last file change: 22-Okt-2012 16:27:00
%
% This program is free software; you can redistribute it and/or modify,
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General gaitPublic License for more details.
%
% You should have received a copy of the GNU General Public License along with this program;
% if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110, USA.
%
% You will find further information about Gait-CAD in the manual or in the following conference paper:
%
% STEGMAIER,J.;ALSHUT,R.;REISCHL,M.;MIKUT,R.: Information Fusion of Image Analysis, Video Object Tracking, and Data Mining of Biological Images using the Open Source MATLAB Toolbox Gait-CAD.
% In:  Proc., DGBMT-Workshop Biosignal processing, Jena, 2012, pp. 109-111; 2012
% Online available: http://www.degruyter.com/view/j/bmte.2012.57.issue-s1-B/bmt-2012-4073/bmt-2012-4073.xml
%
% Please refer to this paper, if you use Gait-CAD with the ImVid extension for your scientific work.

%% check if the tracklets contain additional features, otherwise initialize empty feature list
if (isfield(tracklets, 'featureNames'))
    featureNames = tracklets(1).featureNames;
    numFeatures = size(tracklets(1).features,2);
else
    %% add selected time series of d_orgs to the tracklets
    if (exist('parameter', 'var') && isfield(parameter, 'gui') && isfield(parameter.gui, 'merkmale_und_klassen') && isfield(parameter.gui.merkmale_und_klassen, 'ind_zr'))
        featureIndices = parameter.gui.merkmale_und_klassen.ind_zr;
        featureNames = var_bez(featureIndices, :);
        numFeatures = length(featureIndices);

        for i=1:length(tracklets)
           tracklets(i).features = squeeze(d_orgs(tracklets(i).ids(1), tracklets(i).startTime:tracklets(i).endTime, featureIndices));
           tracklets(i).featureNames = featureNames;
        end
    else
       featureNames = char();
       numFeatures = 0;
    end
end

maxEndTime = 0;
for i=1:length(tracklets)
    if (tracklets(i).endTime > maxEndTime)
        maxEndTime = tracklets(i).endTime;
    end
end

% d_orgs_new = zeros(size(d_orgs,1), size(d_orgs,2), 7+numFeatures);
d_orgs_new = NaN(size(tracklets,2), maxEndTime, 7+numFeatures);

% initialize freeRow variable
freeRow = 1;
% initialize max_free_Row
max_free_Row = 1;

%% iterate over all tracklets and add them to d_orgs
for i=1:length(tracklets)
    
    %% get the current tracklet
    currentTracklet = tracklets(i);
    
    %% identify a free row where the tracklet can be placed in
    if (exist('useCompactDOrgsRepresentation', 'var') && useCompactDOrgsRepresentation == true)
         freeIndices = find(sum(d_orgs_new(:,currentTracklet.startTime:currentTracklet.endTime,1),2) == 0);
         if (isempty(freeIndices))
             freeRow = size(d_orgs_new,1)+1;
         else
             freeRow = min(freeIndices);
         end
    else
        freeRow = i;
    end
    
    % Log the maximal counter for free Row
    if(freeRow > max_free_Row)
        max_free_Row = freeRow;
    end;
        
    %% update the tracklet ids variable and add the tracklet to d_orgs
    tracklets(i).ids = freeRow * ones((currentTracklet.endTime - currentTracklet.startTime + 1), 1)';
    d_orgs_new(freeRow, currentTracklet.startTime:currentTracklet.endTime, 1) = i;
    d_orgs_new(freeRow, currentTracklet.startTime:currentTracklet.endTime, 2) = currentTracklet.endTime - currentTracklet.startTime + 1;
    d_orgs_new(freeRow, currentTracklet.startTime:currentTracklet.endTime, 3:5) = currentTracklet.pos;
    
    if (numFeatures ~= 0)
        d_orgs_new(freeRow, currentTracklet.startTime:currentTracklet.endTime, 6:(6+numFeatures-1)) = currentTracklet.features;
    end
    
    featureIndexOffset = 0;
    if (isfield(currentTracklet, 'state'))
        featureIndexOffset = 1;
        d_orgs_new(freeRow, currentTracklet.startTime:currentTracklet.endTime, (6+numFeatures)) = currentTracklet.state;
    end
    d_orgs_new(freeRow, currentTracklet.startTime:currentTracklet.endTime, (6+numFeatures+featureIndexOffset)) = 1;
    d_orgs_new(freeRow, currentTracklet.endTime, (6+numFeatures+featureIndexOffset)) = 0;
    d_orgs_new(freeRow, currentTracklet.startTime:currentTracklet.endTime, (6+numFeatures+1+featureIndexOffset)) = i;
    
    if (isfield(tracklets, 'divisionOcurred'))
        d_orgs_new(freeRow, currentTracklet.startTime:currentTracklet.endTime, (6+numFeatures+2+featureIndexOffset)) = currentTracklet.divisionOcurred;
    end
    
    if (mod(i, 100) == 0)
        disp(['Processed ' num2str(i) ' / ' num2str(length(tracklets))]);
    end
end

%% Truncate the empty rows of d_orgs containing no tracklets
d_orgs_new = d_orgs_new(1:max_free_Row,:,:);


%% update d_orgs
d_orgs = d_orgs_new;

if (numFeatures ~= 0)
    if (isfield(tracklets, 'state'))
        var_bez = char('id', 'length', 'xpos', 'ypos', 'zpos', featureNames, 'state', 'Tracking state', 'trackletID');
    else
        var_bez = char('id', 'length', 'xpos', 'ypos', 'zpos', featureNames, 'Tracking state', 'trackletID');
    end
else
    if (isfield(tracklets, 'state'))
        var_bez = char('id', 'length', 'xpos', 'ypos', 'zpos', 'state', 'Tracking state', 'trackletID');
    else
        var_bez = char('id', 'length', 'xpos', 'ypos', 'zpos', 'Tracking state', 'trackletID');
    end
end

%% add time series name for division events
if (isfield(tracklets, 'divisionOcurred'))
    var_bez = char(var_bez, 'divisionOcurred');
end
code = ones(size(d_orgs,1),1);

%% save the d_orgs project
save([parameter.projekt.pfad '\' parameter.projekt.datei '_BackTracked.prjz'], '-v7.3', '-mat', 'd_orgs', 'code', 'var_bez');
save([parameter.projekt.pfad '\' parameter.projekt.datei '_BackTracked.tracklets'], '-v7', '-mat', 'tracklets', 'trackletsPerTimePoint');
disp(['Finished generating d_orgs from tracklets. Final project saved to: ' parameter.projekt.pfad '\' parameter.projekt.datei '_BackTracked.prjz']);
%
% next_function_parameter = [parameter.projekt.pfad '\' parameter.projekt.datei '_BackTracked.prjz'];
% ldprj_g;
% callback_load_tracklets;