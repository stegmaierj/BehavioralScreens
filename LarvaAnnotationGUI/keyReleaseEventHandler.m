%%
 % CellCountingGUI.
 % Copyright (C) 2019 J. Stegmaier
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
 
%% the key event handler
function keyReleaseEventHandler(~,evt)
    global settings;

    settings.xLim = get(gca, 'XLim');
    settings.yLim = get(gca, 'YLim');

    %% switch between the images of the loaded series
    if (strcmp(evt.Key, 'rightarrow'))
        settings.currentImage = min(settings.currentImage+1, settings.numImages);
        updateVisualization;
    elseif (strcmp(evt.Key, 'leftarrow'))
        settings.currentImage = max(settings.currentImage-1, 1);
        updateVisualization;
        %% not implemented yet, maybe use for contrast or scrolling
    elseif (strcmp(evt.Character, '+') || strcmp(evt.Key, 'uparrow'))
        settings.gamma(1) = min(5, settings.gamma(1)+0.1);
        updateVisualization;
    elseif (strcmp(evt.Character, '-') || strcmp(evt.Key, 'downarrow'))
        settings.gamma(1) = max(0, settings.gamma(1) - 0.1);
        updateVisualization;
        
    %% toggle correct aspect ratio
    elseif (strcmp(evt.Character, 'a'))
        settings.axesEqual = ~settings.axesEqual;
        updateVisualization;
            
    %% select a region of interest and assign the current group to contained detections
    elseif (strcmp(evt.Character, 's'))
        
      SaveProject();
        
    %% reset the zoom
    elseif (strcmp(evt.Character, 'o'))
        settings.xLim = [1, size(settings.rawImage{settings.currentImage},1)];
        settings.yLim = [1, size(settings.rawImage{settings.currentImage},2)];
        updateVisualization;

    %% show the help dialog
    elseif (strcmp(evt.Character, 'h'))    
        showHelp;
    end
end