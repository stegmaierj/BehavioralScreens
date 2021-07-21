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
 
function mouseUp(obj, evt)

    %% add the settings variable
    global settings;
    
    modifiers = get(gcf,'currentModifier');
    shiftPressed = ismember('shift',modifiers);

    currentPoint = get(gca, 'CurrentPoint');
    currentPoint = round(currentPoint(1,1:2));
            
    %% handle different mouse up events
    selectionType = get(obj, 'SelectionType');
    
    
    %% left pressed
    if (strcmp(selectionType, 'normal') || strcmp(selectionType, 'extend'))
        
        if (isempty(settings.currentDetections{settings.currentImage}))
           settings.markTail = false; 
        end

        if (settings.markTail == true)
            settings.currentDetections{settings.currentImage}(end, settings.tailPosIndices) = currentPoint;
            settings.markTail = false;
            disp(['Adding tail detection at ' num2str(currentPoint)]);

            currentHeadPos = round(settings.currentDetections{settings.currentImage}(end, settings.headPosIndices));
            currentTailPos = round(settings.currentDetections{settings.currentImage}(end, settings.tailPosIndices));

            %% perform livewire tracing
            if (~shiftPressed)
                maxDistance = 50;

                rangeX = (min(currentHeadPos(2), currentTailPos(2))):max(currentHeadPos(2), currentTailPos(2));
                rangeY = (min(currentHeadPos(1), currentTailPos(1))):max(currentHeadPos(1), currentTailPos(1));
                currentPatch = imadjust(settings.rawImage{settings.currentImage}(rangeX, rangeY));

                costImage = fLiveWireGetCostFcn(255 * (max(currentPatch(:)) - currentPatch));
                [iPX, iPY] = fLiveWireCalcP(costImage, currentHeadPos(1) - rangeY(1) + 1, currentHeadPos(2) - rangeX(1) + 1, maxDistance);


                 %% get the live wire prediction
                [iXPath, iYPath] = fLiveWireGetPath(iPX, iPY, currentTailPos(1) - rangeY(1) + 1, currentTailPos(2) - rangeX(1) + 1);
                iXPath = iXPath + rangeY(1) - 1;
                iYPath = iYPath + rangeX(1) - 1;
            else
                [iXPath, iYPath] = bresenham(currentHeadPos(1), currentHeadPos(2), currentTailPos(1), currentTailPos(2));
            end                                

            numPoints = length(iXPath);
            settings.currentLarva{settings.currentImage}(end, 1:numPoints, 1) = iXPath;
            settings.currentLarva{settings.currentImage}(end, 1:numPoints, 2) = iYPath;

        else
            %% update the current detections
            settings.currentDetections{settings.currentImage} = [settings.currentDetections{settings.currentImage}; [size(settings.currentDetections{settings.currentImage},1), 1, currentPoint(1), currentPoint(2), 1, currentPoint(1), currentPoint(2), 1]];

            %% refine centroid for head region
            headRadius = 3;
            imageSize = size(settings.rawImage{settings.currentImage});
            rangeX = max(1, (currentPoint(2)-headRadius)):min(imageSize(1), (currentPoint(2)+headRadius));
            rangeY = max(1, (currentPoint(1)-headRadius)):min(imageSize(2), (currentPoint(1)+headRadius));
            currentHeadPatch = settings.rawImage{settings.currentImage}(rangeX, rangeY);
            [headX, headY] = ind2sub(size(currentHeadPatch), find(currentHeadPatch <= graythresh(currentHeadPatch(:))));

            refinedHeadPosition = currentPoint + ([mean(headY)-headRadius-1, mean(headX)-headRadius-1]);
            settings.currentDetections{settings.currentImage}(end, settings.headPosIndices) = refinedHeadPosition;

            if (isempty(settings.currentLarva{settings.currentImage}))
                settings.currentLarva{settings.currentImage} = zeros(1, 0, 2);
            else
                settings.currentLarva{settings.currentImage}(size(settings.currentDetections{settings.currentImage},1), :, :) = 0;
            end

            settings.markTail = true;
            disp(['Adding head detection at ' num2str(currentPoint)]);
        end
            
           
        %% delete detection closest to the click position
    elseif (strcmp(selectionType, 'alt'))
            
        %% get the current detections and the current click position
        currentDetections = settings.currentDetections{settings.currentImage};

        %% identify the closest existing detection to the click position
        distances = sqrt((currentDetections(:,settings.headPosIndices(1)) - currentPoint(1)).^2 + (currentDetections(:,settings.headPosIndices(2)) - currentPoint(2)).^2);
        [~, minIndex] = min(distances);

        %% update the current detections
        settings.currentDetections{settings.currentImage}(minIndex, :) = [];
        settings.currentLarva{settings.currentImage}(minIndex, :, :) = [];
        disp(['Removing detection at ' num2str(currentPoint)]);
                 
    end    
    
    %% update the visualization
    updateVisualization;
end