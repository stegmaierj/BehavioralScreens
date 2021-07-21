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
 
inputFolder = 'C:\Users\stegmaier\Downloads\Labeling\RawVideos\';
inputFiles = dir([inputFolder '*.avi']);
outputFolder = 'C:\Users\stegmaier\Downloads\Labeling\RawImages\';
numImagesPerSequence = 10;

for f=1:length(inputFiles)

    inputFile = [inputFolder inputFiles(f).name];
    inputVideo = VideoReader( inputFile );

    numFrames = round(inputVideo.Duration * inputVideo.FrameRate);
    
    stepSize = round(numFrames / numImagesPerSequence);
    
    frameRange = 1:stepSize:numFrames;
    
    
    numFramesAnalysis = length(frameRange);

    for i=frameRange
        currentFrame = im2double(read(inputVideo, i));
        if (size(currentFrame,3) > 2)
            currentFrame = rgb2gray(currentFrame);
        end

        imwrite(currentFrame, sprintf('%s%s_Frame=%05d.png', outputFolder, strrep(inputFiles(f).name, '.avi', ''), i));
    end
end