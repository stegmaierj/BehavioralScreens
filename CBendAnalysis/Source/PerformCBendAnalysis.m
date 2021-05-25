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

function [] = PerformCBendAnalysis(inputFile, outputDir, settings)

    %% setup required parameters
    if (~isfield(settings,'roiRadiusPadding')); roiRadiusPadding = 0; else; roiRadiusPadding = settings.roiRadiusPadding; end %% used to further reduce the circular roi, i.e., the roi mask radius is shrinked by this value
    if (~isfield(settings,'minArea')); minArea = 10; else; minArea = settings.minArea; end %% used to filter out small segments such as background noise
    if (~isfield(settings,'maxArea')); maxArea = 100; else; maxArea = settings.maxArea; end %% used to filter large segments such as parts of the border regions
    if (~isfield(settings,'frameRate')); frameRate = 1000; else; frameRate = settings.frameRate; end %% the number of frames per second (used to scale latency values to ms)
    if (~isfield(settings,'minPeakHeight')); minPeakHeight = 20; else; minPeakHeight = settings.minPeakHeight; end %% minimium angle change in degree. Only c-bends with head angle changes larger than this value are considered a valid c-bend
    if (~isfield(settings,'minPeakWidth')); minPeakWidth = 3; else; minPeakHeight = settings.minPeakHeight; end %% minimium angle change in degree. Only c-bends with head angle changes larger than this value are considered a valid c-bend
    if (~isfield(settings,'minPeakProminence')); minPeakProminence = 20; else; minPeakProminence = settings.minPeakProminence; end %% minimum relative height of a peak, can be used to supress small peaks between two larger bends
    if (~isfield(settings,'minTrackletLengthAfterStimulus')); minTrackletLengthAfterStimulus = 50; else; minTrackletLengthAfterStimulus = settings.minTrackletLengthAfterStimulus; end %% tracklets that are shorter than this minimum length are discarded, as they potentially do not contain both c-bends
    if (~isfield(settings,'framePadding')); framePadding = 50; else; framePadding = settings.framePadding; end %% this amount of frames is used before and after the led-onset, i.e. it should be large enough to cover the desired c-bend behaviors
    %if (~isfield(settings,'larvaLength')); larvaLength = 45; else; larvaLength = settings.larvaLength; end % length of a straight larva in pixels (currently not used)
    if (~isfield(settings,'startFrame')); startFrame = 0; else; startFrame = settings.startFrame; end  %% skips all frames before that frame to reduce the processing time
    if (~isfield(settings,'numFramesAfterStimulus')); numFramesAfterStimulus = 1000; else; numFramesAfterStimulus = settings.numFramesAfterStimulus; end  %% the number of frames to consider after the pulse
    if (~isfield(settings,'generateResultVideo')); generateResultVideo = false; else; generateResultVideo = settings.generateResultVideo; end %% if enabled, result videos are generated
    if (~isfield(settings,'latencyThreshold')); latencyThreshold = 16; else; latencyThreshold = settings.latencyThreshold; end %% head orientation changes above this value are considered as cbend start (16° as described in Burgess07)
    if (~isfield(settings,'minMaxRatioLEDIntensity')); minMaxRatioLEDIntensity = 1.01; else; minMaxRatioLEDIntensity = settings.minMaxRatioLEDIntensity; end %% minimum vs. maximum intensity ratio of the led. If close to 1 means that the LED-related intensity increase is hardly visible.
    
    if (~isfield(settings,'xRangeLED')); xRangeLED = 1:100; else; xRangeLED = settings.xRangeLED; end %% range of pixels in the x-direction to locate the LED
    if (~isfield(settings,'yRangeLED')); yRangeLED = 1:100; else;yRangeLED = settings.yRangeLED; end %% range of pixels in the y-direction to locate the LED
    
    if (~isfield(settings,'debugFigures')); debugFigures = false; else; debugFigures = settings.debugFigures; end %% if enabled, debug figures are shown

    %% specify the input file and output directory
    [~, file, ~] = fileparts(inputFile);
    logFile = fopen([outputDir '\' file '_Logfile.txt'], 'wt');

    %% load the input video
    fboth(logFile, '- Trying to open input video ...'); tic;
    inputVideo = VideoReader( inputFile );

    if (~exist(inputFile, 'file'))
        fboth(logFile, ['ERROR: The specified input file ' strrep(inputFile, '\', '\\') ' does not exist or cannot be read! Please specify a different file or install an appropriate codec.']);
    else
        fboth(logFile, ['- Successfully opened input video ' strrep(inputFile, '\', '\\') ' for reading ...']);
    end

    %% extract the frames and write the raw images to disk
    loadInputVideo = true;
    if (loadInputVideo == true)
    fboth(logFile, '+ Calculating mean image and detecting pulse frame ranges ...'); tic;
    numFrames = round(inputVideo.Duration * inputVideo.FrameRate)-startFrame;
    meanImage = zeros(size(inputVideo.Height, inputVideo.Width));
    meanIntensityLED = zeros(numFrames,1);
    frameRange = 1:numFrames;
    numFramesAnalysis = length(frameRange);

    for i=frameRange
        currentFrame = im2double(read(inputVideo, i+startFrame));
        if (size(currentFrame,3) > 2)
            currentFrame = rgb2gray(currentFrame);
        end

        meanImage = meanImage + (1/numFramesAnalysis) * currentFrame;
        meanIntensityLED(i) = mean(mean(currentFrame(yRangeLED, xRangeLED)));
        if (mod(i, 100) == 0 || i == numFrames)
            fboth(logFile, ['\t- Processed ' num2str(round(100*i/numFrames)) '%%']);
        end
    end

    %% identify the led on/off time points
    binarizedMeanIntensityLED = meanIntensityLED > (0.5*(quantile(meanIntensityLED, 0.05)+quantile(meanIntensityLED, 0.95)));
    ledSwitchTimePoints = find(binarizedMeanIntensityLED - circshift(binarizedMeanIntensityLED, 1));
    numPulses = floor(length(ledSwitchTimePoints)/2);
    
    minMaxRatio = quantile(meanIntensityLED, 0.95) / quantile(meanIntensityLED, 0.05);
    
    if (minMaxRatio < minMaxRatioLEDIntensity)
        fboth(logFile, ['- No pulses were detected in the image. If there is indeed a pulse present, try to decrease the minMaxRatioLEDIntensity parameter (currently set to ' num2str(minMaxRatioLEDIntensity) ') or relocate the ROI for the LED using the parameters xRangeLED and yRangeLED.']);
        numPulses = 0;
        ledSwitchTimePoints = [];
    end
    
    frameRanges = zeros(numPulses, 2);
    for i=1:2:length(ledSwitchTimePoints)
        frameRanges(((i+1)/2),:) = ledSwitchTimePoints(i:i+1)';
        fboth(logFile, sprintf('- Pulse %i ranges from %i to %i', (i+1)/2, ledSwitchTimePoints(i), ledSwitchTimePoints(i+1)));
    end
    fboth(logFile, ['- Finished calculating mean image and pulse detection in ' num2str(toc) ' seconds']);

    %% extract the well mask
    circleROIImage = ExtractWellMask(meanImage, roiRadiusPadding);
    end

    tempImages = struct();

    %% separately analyze the different pulse sequences
    for p=1:numPulses

        %% initialize larva models and frame range etc.
        larvaModels = cell(numFrames, 1);
        %skeletonResults = cell(numFrames, 1);
        %frameRange = (frameRanges(p,1)-framePadding):(frameRanges(p,1)+framePadding);

        frameRange = max(1, (frameRanges(p,1)-framePadding)):min(numFrames, (frameRanges(p,1)+numFramesAfterStimulus));

        parameter.projekt.pfad = outputDir;
        parameter.projekt.datei = [file '_Pulse' num2str(p)];

        fboth(logFile, sprintf('------------------- Processing Pulse %i (frames:%i-%i) -------------------', p, frameRange(1), frameRange(end)));

        fboth(logFile, '+ Performing image pre-processing for larva detection ...'); tic;
        for i=frameRange

            %% get the current frame
            currentFrame = im2double(read(inputVideo, i+startFrame));
            if (size(currentFrame, 3) > 2)
                currentFrame = rgb2gray(currentFrame);
            end

            %% process the input image
            inputImage = imcomplement(imadjust(currentFrame));
            thresholdedImage = adaptivethreshold(inputImage, 30, -0.5*std(inputImage(:))) .* circleROIImage;
            headImage = imopen(thresholdedImage, strel('disk', 2));

            %% calculate skeleton image and discard small objects
            %skeletonImage = bwmorph(thresholdedImage > 0, 'skel', inf);
            tempImages(i).currentLabelImage = bwlabel(headImage | ~circleROIImage);
            currentRegionProps = regionprops(tempImages(i).currentLabelImage, 'Area', 'PixelIdxList');
            for j=1:length(currentRegionProps)
                if (currentRegionProps(j).Area < minArea || currentRegionProps(j).Area > maxArea)
                    tempImages(i).currentLabelImage(currentRegionProps(j).PixelIdxList) = 0;
                end
            end

            %% extract the larva models
            larvaModels{i} = IdentifyLarvaModels(bwlabeln(tempImages(i).currentLabelImage > 0), inputImage, ones(size(circleROIImage)));

            %% show progress (TODO: potentially use a log file for the progress)
            if (mod(i-frameRange(1), 10) == 0 || i == numFrames)
                fboth(logFile, ['\t- Processed ' num2str(round(100*(i-frameRange(1))/(frameRange(end)-frameRange(1)))) '%%']);
            end
        end
        fboth(logFile, ['- Finished larva detection in ' num2str(toc) ' seconds']);

        %% perform tracking of the detected larvae
        fboth(logFile, '+ Performing tracking of all detected larvae ...'); tic;
        d_orgs = zeros(0, length(frameRange), 14);
        var_bez = char('id', 'area', 'xpos', 'ypos', 'zpos', 'headDirX', 'headDirY', 'eccentricity', 'solidity', 'convexArea', 'equivDiameter', 'angleChange', 'relativeAngleSum', 'angleChangeSum');
        numDetectedLarva = 0;
        for i=1:length(frameRange)
            currentFrame = frameRange(i);
            
            %% continue if there's no larva contained in the current frame
            if (~isfield(larvaModels{currentFrame}, 'larvaModel'))
                continue;
            end

            for j=1:length(larvaModels{currentFrame})
               
                %% identify the head direction
                if (size(larvaModels{currentFrame}(j).larvaModel, 1) < 2)
                    headDirX = 0;
                    headDirY = 0;
                else
                    headDirX = larvaModels{currentFrame}(j).larvaModel(2,1) - larvaModels{currentFrame}(j).larvaModel(1,1);
                    headDirY = larvaModels{currentFrame}(j).larvaModel(2,2) - larvaModels{currentFrame}(j).larvaModel(1,2);
                end

                %% save features to d_orgs (only needed for the tracking)
                d_orgs(j,i,:) = [j, ...
                    larvaModels{currentFrame}(j).area, ...
                    larvaModels{currentFrame}(j).larvaModel(1,1), ...
                    larvaModels{currentFrame}(j).larvaModel(1,2), ...
                    0, ...
                    headDirX, ...
                    headDirY, ...
                    larvaModels{currentFrame}(j).eccentricity, ...
                    larvaModels{currentFrame}(j).solidity, ...
                    larvaModels{currentFrame}(j).convexArea, ...
                    larvaModels{currentFrame}(j).equivDiameter, ...
                    sum(abs(larvaModels{currentFrame}(j).angleChange)), ...,
                    larvaModels{currentFrame}(j).relativeAngleSum, ...,
                    larvaModels{currentFrame}(j).angleChangeSum];
                
                %% count total number of larva
                numDetectedLarva = numDetectedLarva + 1;
            end
        end
        
        if (numDetectedLarva == 0)
            fboth(logFile, ['- No larva detected during current pulse, skipping tracking ...']);
            continue;
        end

        %% perform the tracking using a nearest neighbor tracking approach
        parameter.gui.merkmale_und_klassen.ind_zr = [2,6:14];
        code = ones(size(d_orgs,1), 1);
        callback_tracking_regionprops;
        callback_extract_tracklets_regionprops;
        parameter.gui.merkmale_und_klassen.ind_zr = [2,6:12];
        convert_tracklets_structure_to_dorgs;
        fboth(logFile, ['- Finished larva tracking in ' num2str(toc) ' seconds']);

        %% iterate through the tracklets and calculate the cbend latency, angular velocity etc.
        fboth(logFile, '+ Analyzing larva behavior ...'); tic;
        completeTracks = [];
        for j=1:length(tracklets)

            %% convert the relative time range of the tracklets to the absolute time frame
            tracklets(j).startTime = tracklets(j).startTime-framePadding+frameRanges(p,1); %#ok<AGROW>
            tracklets(j).endTime = tracklets(j).endTime-framePadding+frameRanges(p,1); %#ok<AGROW>
            if (tracklets(j).startTime == tracklets(j).endTime || ...
                    tracklets(j).endTime <= frameRanges(p,1) || ...
                    tracklets(j).startTime > frameRanges(p,1) || ...
                    (tracklets(j).endTime - frameRanges(p,1)) < minTrackletLengthAfterStimulus)

                fboth(logFile, ['\t- Skipped tracklet with length ' num2str(length(tracklets(j).ids)) '']);
                continue;
            else
                fboth(logFile, ['\t- Processing tracklet with length ' num2str(length(tracklets(j).ids)) '']);
            end

            %% calculate the angular changes to the start head angle (idle position)
            relativeAngles = [];
            dirXFeatureIndex = 2;
            dirYFeatureIndex = 3;
            %eccentricityIndex = 4;
            relativeAngleIndex = 9;
            startDirection = tracklets(j).features(framePadding,[dirXFeatureIndex, dirYFeatureIndex]);
            startDirection = startDirection / norm(startDirection);

            %% calculate the relative head angle for all frames after the onset of the led/pulse
            relativePulseTimePoint = find(tracklets(j).startTime:tracklets(j).endTime == frameRanges(p,1));
            for k=relativePulseTimePoint:length(tracklets(j).ids)
                currentDirection = tracklets(j).features(k,[dirXFeatureIndex, dirYFeatureIndex]);
                currentDirection = currentDirection / norm(currentDirection);

                %% add angle to the array of all angles
                relativeAngles = [relativeAngles; rad2deg(real(acos(dot(startDirection, currentDirection))))]; %#ok<AGROW>
            end
            tracklets(j).features(relativePulseTimePoint:length(tracklets(j).ids),relativeAngleIndex) = relativeAngles; %#ok<AGROW>
            tracklets(j).featureNames = char(tracklets(j).featureNames, 'relativeAngles'); %#ok<AGROW>
            tracklets(j).relativePulseTimePoint = relativePulseTimePoint; %#ok<AGROW>
            tracklets(j).impulseStartTimePoint = frameRanges(p,1); %#ok<AGROW>
            tracklets(j).impulseEndTimePoint = frameRanges(p,2); %#ok<AGROW>

            frameToFrameDistance = sqrt(sum((circshift(tracklets(j).pos(:,1:2), -1) - tracklets(j).pos(:,1:2)).^2, 2));
            tracklets(j).distanceTraveledAfterPulse = sum(frameToFrameDistance(1:(end-1))); %#ok<AGROW>

            %% extract the peaks
            [peakValues, peakIndices] = findpeaks(smooth(relativeAngles, 1, 'moving'), 'MinPeakHeight', minPeakHeight, 'MinPeakWidth', minPeakWidth, 'MinPeakProminence', minPeakProminence);

            %% identify the latency, starting from the first peak, search for the position where the angular change drops first below the threshold value
            belowThresholdAngles = find(relativeAngles < latencyThreshold);

            %% save the first two peaks as the first two cbend events and calculate latency, angles and angular velocity
            if (length(peakIndices) >= 1)
                currentLatencyIndex = max(belowThresholdAngles(find(belowThresholdAngles < peakIndices(1)))); %#ok<FNDSB>
                currentLatency = 1/frameRate * currentLatencyIndex * 1000;

                %% calculate the maximum angular velocity in °/ms
                angleChanges1 = circshift(relativeAngles(currentLatencyIndex:peakIndices(1)), -1) - relativeAngles(currentLatencyIndex:peakIndices(1));
                maxAngularVelocity1 = max(angleChanges1(1:end-1)) / (1/frameRate) / 1000;

                tracklets(j).latency = currentLatency; %#ok<AGROW>
                tracklets(j).firstPeakDistance = peakIndices(1); %#ok<AGROW>
                tracklets(j).firstPeakTime = 1/frameRate * tracklets(j).firstPeakDistance * 1000; %#ok<AGROW>
                tracklets(j).firstPeakDuration = tracklets(j).firstPeakTime - currentLatency; %#ok<AGROW>
                tracklets(j).firstPeakAngle = peakValues(1); %#ok<AGROW>
                tracklets(j).firstPeakAngleMaxAngularVelocity = maxAngularVelocity1; %#ok<AGROW>
            end

            if (length(peakIndices) >= 2)
                [~, minIndexBetweenPeaks] = min(relativeAngles(peakIndices(1):peakIndices(2)));

                %% calculate the maximum angular velocity in °/ms
                angleChanges2 = circshift(relativeAngles((peakIndices(1)+minIndexBetweenPeaks):peakIndices(2)), -1) - relativeAngles((peakIndices(1)+minIndexBetweenPeaks):peakIndices(2));
                maxAngularVelocity2 = max(angleChanges2(1:end-1)) / (1/frameRate) / 1000;

                tracklets(j).secondPeakDistance = peakIndices(2); %#ok<AGROW>
                tracklets(j).secondPeakTime = 1/frameRate * tracklets(j).secondPeakDistance * 1000; %#ok<AGROW>
                tracklets(j).secondPeakDuration = tracklets(j).secondPeakTime - (peakIndices(1) + minIndexBetweenPeaks); %#ok<AGROW>
                tracklets(j).secondPeakAngle = peakValues(2); %#ok<AGROW>
                tracklets(j).secondPeakAngleMaxAngularVelocity = maxAngularVelocity2; %#ok<AGROW>
            end

            if ((framePadding+numFramesAfterStimulus) == size(tracklets(j).pos,1))
                completeTracks = [completeTracks; j]; %#ok<AGROW>
            end

            %%
            if (debugFigures == true)
                figure;
                plot(smooth(relativeAngles, 1, 'moving')); hold on;
                plot(peakIndices, peakValues, '*g');
            end
        end
        
        %% delete intermediate results
        allFilesAndFolders = dir(outputDir);
        justFilesInDir = allFilesAndFolders(~([allFilesAndFolders.isdir]));
        numberOfFiles = length(justFilesInDir);
        for i = 1:numberOfFiles
            if strfind(justFilesInDir(i).name, 'prjz')
                delete([outputDir '/' justFilesInDir(i).name]);
            elseif strfind(justFilesInDir(i).name, 'tracklets')
                delete([outputDir '/' justFilesInDir(i).name]);
            end
        end

        %% save the results of the current pulse
        save([outputDir filesep file '_Pulse' sprintf('%03d', p) '.mat'], 'tracklets', 'larvaModels', 'd_orgs', 'code', 'var_bez', 'completeTracks');

        %% display result message
        fboth(logFile, ['- Finished larva feature analysis in ' num2str(toc) ' seconds and saved results to ' strrep([parameter.projekt.pfad parameter.projekt.datei], '\', '\\') '.mat']);

        %% generate result videos
        if (generateResultVideo == true)
            fh = figure('Visible', 'off'); clf;
            colordef black; %#ok<COLORDEF>
            set(gcf, 'color', 'black');

            writerObj = VideoWriter([outputDir filesep file '_Pulse' num2str(p) '.avi']); %#ok<TNMLP>
            writerObj.FrameRate = 15;
            open(writerObj);

            %% plot the result figure
            %lastEccentricity = 1;
            for i=frameRange
                %% identify the end points
                clf(fh);
                currentFrame = im2double(read(inputVideo, i+startFrame));
                if (size(currentFrame, 3) > 2)
                    currentFrame = rgb2gray(currentFrame);
                end
                inputImage = imcomplement(imadjust(currentFrame));
                hold on;
                imagesc(inputImage); colormap gray;

                [xpos, ypos] = ind2sub(size(tempImages(i).currentLabelImage), find(tempImages(i).currentLabelImage>0));
                plot(ypos, xpos, '.b');
                axis equal;


                for j=1:length(tracklets)
                    currentFrameRelative = (i-tracklets(j).startTime+1);
                    if (size(tracklets(j).pos,1) >= currentFrameRelative && currentFrameRelative > 0)
                        plot(tracklets(j).pos(1:currentFrameRelative,1), tracklets(j).pos(1:currentFrameRelative,2), '-r');
                        text(tracklets(j).pos(currentFrameRelative,1), tracklets(j).pos(currentFrameRelative,2), ['ID = ' num2str(j)]);
                    end
                end
                hold off;

                fboth(logFile, sprintf('- Processed video frame %i of %i', i, frameRange(end)));
                writeVideo(writerObj,getframe(fh));
            end
            close(fh);
            close(writerObj);
        end
    end

    fclose(logFile);
end