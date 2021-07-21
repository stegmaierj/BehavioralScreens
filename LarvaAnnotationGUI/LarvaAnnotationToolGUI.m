%%
 % LarvaAnnotationToolGUI.
 % Copyright (C) 2021 J. Stegmaier
 %
 % Licensed under the Apache License, Version 2.0 (the "License");
 % you may not use this file except in compliance with the License.
 % You may obtain a copy of the Liceense at
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

%% initialize the global settings variable
close all;
if (exist('settings', 'var'))
    clearvars -except settings;
else
    clearvars;
end
global settings;

%% add paths for the live wire and the image reader
if (~isdeployed)
    addpath('ThirdParty/LiveWire/');
    addpath('ThirdParty/bresenham/');
    addpath('ThirdParty/');
end

%% open the input image
settings.inputPath = uigetdir('', 'Select the folder containing the images for counting.');
settings.inputPath = [settings.inputPath filesep];

settings.outputPath = [settings.inputPath '..' filesep];
settings.outputPathLabelImages = [settings.outputPath 'LabelImages' filesep];
settings.outputPathMaskImages = [settings.outputPath 'HeadImages' filesep];
settings.outputPathDontCareImages = [settings.outputPath 'DontCareImages' filesep];
if (~isfolder(settings.outputPathLabelImages)); mkdir(settings.outputPathLabelImages); end
if (~isfolder(settings.outputPathMaskImages)); mkdir(settings.outputPathMaskImages); end
if (~isfolder(settings.outputPathDontCareImages)); mkdir(settings.outputPathDontCareImages); end

%% identify suitable files and initialize the results table
settings.imageFiles = dir([settings.inputPath filesep '*.png']);
settings.numImages = length(settings.imageFiles);
settings.resultsTable = cell(settings.numImages, 2);

%% load  the input images
for i=1:settings.numImages
    settings.rawImagePath{i} = [settings.inputPath settings.imageFiles(i).name];
    settings.rawImage{i} = im2double(imread(settings.rawImagePath{i}));
    if (size(settings.rawImage{i},3) > 2)
        settings.rawImage{i} = rgb2gray(settings.rawImage{i});
    end
    
    settings.resultsTable{i} = settings.rawImagePath{i};
end
settings.currentDetections = cell(settings.numImages, 1);
settings.currentLarva = cell(settings.numImages, 1);

%% check if previous labelings exist and load if existent
if (exist([settings.outputPath 'detections.mat'], 'file'))
    load([settings.outputPath 'detections.mat']);
    settings.currentDetections = currentDetections;
    settings.currentLarva = currentLarva;
end

%% initialize the settings
settings.currentImage = 1;
settings.colormapIndex = 1;
settings.markerSize = 15;
settings.gamma = 1;
settings.minIntensity = min(settings.rawImage{settings.currentImage}(:));
settings.maxIntensity = max(settings.rawImage{settings.currentImage}(:));
settings.axesEqual = false;
settings.fontSize = 14;
settings.colormapStrings = {'gray', 'parula', 'jet'};

settings.headPosIndices = 3:4;
settings.tailPosIndices = 6:7;
settings.markTail = false;

%% specify the figure boundaries
settings.xLim = [0, size(settings.rawImage{settings.currentImage},1)];
settings.yLim = [0, size(settings.rawImage{settings.currentImage},2)];


%% open the main figure
settings.mainFigure = figure(1);

%% mouse, keyboard events and window title
set(settings.mainFigure, 'WindowScrollWheelFcn', @scrollEventHandler);
set(settings.mainFigure, 'KeyReleaseFcn', @keyReleaseEventHandler);
set(settings.mainFigure, 'WindowButtonDownFcn', @mouseUp);
% set(settings.mainFigure, 'WindowButtonMotionFcn', @mouseMove);
set(settings.mainFigure, 'CloseRequestFcn', @closeRequestHandler);

%% update the visualization
updateVisualization;