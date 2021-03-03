%%
% CBendAnalysis.
% Copyright (C) 2021 R. Peravali, D. Marcato, R. Mikut, J. Stegmaier
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Please refer to the documentation for more information about the software
% as well as for installation instructions.
%
% If you use this application for your work, please cite the repository and one
% of the following publications:
%
% TBA
%
%%

%% add dependencies
addpath('ThirdParty/');

%% open input folder and find all contained avi files
inputFolder = uigetdir;
fileFilter = [inputFolder '/*.avi'];
avifiles = dir(fileFilter);

%% set to true if result videos should be generated
generateResultVideo = false;

%% process all detected videos
parfor i=1:size(avifiles,1)
    
    %% check if video was analyzed before and skip in that case (delete result files to reprocess)
    [pathstr,name,ext] = fileparts([inputFolder '/' avifiles(i).name]);
    if exist([pathstr '/' name '.mat'], 'file')
        continue;
    end
    
    %% otherwise analyse video
    AnalyzeSingleVideo([inputFolder filesep avifiles(i).name], generateResultVideo);
end

%% aggregate all final results to a single table
fileID = fopen([inputFolder filesep 'combinedResults.csv'], 'wb');
fprintf(fileID, 'Filename;Pulse;NumTrackedLarva;NumActiveLarva;MeanLatency;MeanDistanceTraveled;MeanDistanceTraveledActiveLarva\n');

%% find all valid result files
resultFiles = dir([inputFolder '/*Pulse*.mat']);
for i=1:length(resultFiles)
    
    %% load the current results
    currentFile = [inputFolder filesep resultFiles(i).name];
    load(currentFile);
    
    %% open result file for the current video/pulse
    fileIDPulse = fopen([inputFolder filesep strrep(resultFiles(i).name, '.mat', '.csv')], 'wb');
    fprintf(fileIDPulse, 'Filename;Pulse;LarvaID;IsActive;Latency;DistanceTraveled;\n');
    
    %% extract the current summary statistics
    numTrackedLarva = length(completeTracks);
    numActiveLarva = 0;
    meanLatency = 0;  
    meanDistanceTraveled = 0;
    meanDistanceTraveledActiveLarva = 0;
    for j=completeTracks'
        
        currentDistanceTraveledAfterPulse = 0;
        currentLatency = 0;
        isActive = 0;
        
        if (isfield(tracklets, 'distanceTraveledAfterPulse') && ~isempty(tracklets(j).distanceTraveledAfterPulse))
            meanDistanceTraveled = meanDistanceTraveled + tracklets(j).distanceTraveledAfterPulse;
        end
        
        if (isfield(tracklets, 'latency') && ~isempty(tracklets(j).latency))
            isActive = 1;
            currentLatency = tracklets(j).latency;
            currentDistanceTraveledAfterPulse = tracklets(j).distanceTraveledAfterPulse;
            numActiveLarva = numActiveLarva + 1;
            meanLatency = meanLatency + currentLatency;
            meanDistanceTraveledActiveLarva = meanDistanceTraveledActiveLarva + currentDistanceTraveledAfterPulse;
        end
        
        %% write results of the current larva
        fprintf(fileIDPulse, '%s;%i;%i;%i;%.2f;%.2f\n', resultFiles(i).name, str2double(resultFiles(i).name(end-6:end-4)), j, isActive, currentLatency, currentDistanceTraveledAfterPulse);
    end
    
    %% close current pulse results
    fclose(fileIDPulse);
    
    %% scale the measures according to the number of (active) larva
    meanLatency = meanLatency / numActiveLarva;
    meanDistanceTraveled = meanDistanceTraveled / numTrackedLarva;
    meanDistanceTraveledActiveLarva = meanDistanceTraveledActiveLarva / numActiveLarva;
    
    %% write the current result line to the file
    fprintf(fileID, '%s;%i;%i;%i;%.2f;%.2f;%.2f\n', resultFiles(i).name, str2double(resultFiles(i).name(end-6:end-4)), numTrackedLarva, numActiveLarva, meanLatency, meanDistanceTraveled, meanDistanceTraveledActiveLarva);
end

%% close the results file
fclose(fileID);

