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

%% get the global settings variable
global settings;

%% iterate over all images and add the file name as well as the current counts
for i=1:settings.numImages
    
    currentMaskImage = zeros(size(settings.rawImage{i}));
    currentHeadImage = zeros(size(settings.rawImage{i}));
    dontCareImage = ones(size(settings.rawImage{i}));
    [folder, file, ext] = fileparts(settings.rawImagePath{i});
    
    if (~isempty(settings.currentLarva{i}))
    
        for l=1:size(settings.currentLarva{i}, 1)
           currentXPos = settings.currentLarva{i}(l,:,1);
           currentYPos = settings.currentLarva{i}(l,:,2);
           currentXPos = currentXPos(currentXPos > 0);
           currentYPos = currentYPos(currentYPos > 0);

           currentHeadImage(settings.currentLarva{i}(l,1,1), settings.currentLarva{i}(l,1,2)) = l;

           for p=1:length(currentXPos)
              currentMaskImage(currentXPos(p), currentYPos(p)) = l;
           end
        end

        currentHeadImage = imdilate(currentHeadImage, strel('sphere', 3))';
        currentMaskImage = imdilate(currentMaskImage, strel('sphere', 1))';
        currentMaskImage = max(currentMaskImage, currentHeadImage);
        dontCareImage = imdilate(currentMaskImage > 0, strel('sphere', 3)) - (currentMaskImage > 0);
        dontCareImage = ~dontCareImage;
    end
    
    imwrite(uint8(currentMaskImage), [settings.outputPath 'LabelImages' filesep file '_LabelImage.png']);
    imwrite(uint8(currentHeadImage), [settings.outputPath 'HeadImages' filesep file '_HeadImage.png']);
    imwrite(uint8(dontCareImage), [settings.outputPath 'DontCareImages' filesep file '_DontCareImage.png']);
end

%% save the current detections to restore the project
currentDetections = settings.currentDetections;
currentLarva = settings.currentLarva;
save([settings.outputPath 'detections.mat'], 'currentDetections', 'currentLarva');
disp('Successfully saved the project!');