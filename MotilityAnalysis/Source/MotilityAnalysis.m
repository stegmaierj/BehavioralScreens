%% function to analyze a video file containing moving zebrafish larvae
%% inputVideoFileName: the file name of the video to process
%% settings: struct containing the optional parameters. Otherwise default parameters are used

%% RESULTS:
%% averageResults.activePercentage: returns the percentage of movements higher than the "movementThreshold" parameter
%% averageResults.idlePercentage: returns 100 - activePercentage
%% averageResults.totalDistancePixel: total traveled distance in pixels normalized by the number of objects, i.e., by settings.numObjects
%% averageResults.totalDistancePixelMM: total distance traveled in mm --> Be sure to set the settings.pixelToMMFactor parameter correctly for correct scaling.
%% averageResults.swimSpeedPixelPerSecond: swim speed in pixel per second
%% averageResults.swimSpeedMMPerSecond: swim speed in mm per second

function [averageResults, resultImage, d_orgs] = MotilityAnalysis(inputVideoFileName, settings)

%% setup default parameters
if (~isfield(settings, 'minArea')); settings.minArea = 1000; end; %% the minimum area for an object to be considered as valid. Mainly used to supress noise detections
if (~isfield(settings, 'objectRadius')); settings.objectRadius = 50; end; %% potentially splits touching objects if the radius is closer than 50 (experimental...)
if (~isfield(settings, 'numObjects')); settings.numObjects = 1; end; %% the number of objects, used to limit the objects reported in the result tables
if (~isfield(settings, 'numFrames')); settings.numFrames = -1; end; %% number of frames to process, -1 processes the entire video
if (~isfield(settings, 'scaleFactor')); settings.scaleFactor = 1.0; end; %% default is 1.0, if smaller, the processing is performed on a resized version of the video
if (~isfield(settings, 'debugFigures')); settings.debugFigures = false; end; %% enables debug figure/video generation
if (~isfield(settings, 'writeAvi')); settings.writeAvi = false; end; %% if enabled a result movie is created. Note that resultAviFileName is used for the output
if (~isfield(settings, 'resultAviFileName')); settings.resultAviFileName = 'resultVideo.mp4'; end; %% file name for debug output videos
if (~isfield(settings, 'maxAreaRatio')); settings.maxAreaRatio = 1.8; end;  %% if the ratio of the areas of two linked objects is larger than this one, the track is not considered
if (~isfield(settings, 'maxVelocity')); settings.maxVelocity = 30; end; %% maximally allowed distance of two nearest neighbors in pixels (to prevent large jumps)
if (~isfield(settings, 'pixelToMMFactor')); settings.pixelToMMFactor =  10 * 3.4 / 417; end; %% conversion factor from pixel to mm. E.g. measured 417 pixel at the innermost circle, which corresponds to 3.4cm in reality, i.e. 34mm
if (~isfield(settings, 'gaussianSigma')); settings.gaussianSigma = (sqrt(15)*settings.scaleFactor)^2; end; %% gaussian smoothing of the difference image for improved segmentation boundaries
if (~isfield(settings, 'startFrame')); settings.startFrame = 1; end; %% start frame, e.g. if only a subset of the video is of importance
if (~isfield(settings, 'movementThreshold')); settings.movementThreshold = 3; end;   %% smaller movements than this value are not counted for the average movement
if (~isfield(settings, 'minimumIntensityThreshold')); settings.minimumIntensityThreshold = 0.05; end; %% minimum threshold used for the binarization of the difference image
if (~isfield(settings, 'invertIntensity')); settings.invertIntensity = false; end; %% minimum threshold used for the binarization of the difference image
if (~isfield(settings, 'usePartialMeanImage')); settings.usePartialMeanImage = false; end; %% if set to true, only a mean image surrounding the current frame in a radius of "partialMeanImageRadius" is used. Otherwise, the entire sequence is used.
if (~isfield(settings, 'partialMeanImageRadius')); settings.partialMeanImageRadius = 10; end; %% radius for the partial mean image calculation. Only used if "usePartialMeanImage" is set to true.
if (~isfield(settings, 'useLocalAdaptiveThreshold')); settings.useLocalAdaptiveThreshold = true; end; %% uses local adaptive threshold with window size "localAdaptiveThresholdWindowSize". If false, global threshold with Otsu and minimumIntensityThreshold as lower boundary is used.
if (~isfield(settings, 'localAdaptiveThresholdWindowSize')); settings.localAdaptiveThresholdWindowSize = 50; end; %% window size for the local adaptive threshold.
if (~isfield(settings, 'useWatershedBasedSplitting')); settings.useWatershedBasedSplitting = false; end; %% window size for the local adaptive threshold.

%% specify filter mask for gaussian smoothing
filterMask = fspecial('gaussian', 4*[settings.gaussianSigma, settings.gaussianSigma], settings.gaussianSigma);

%% start time measurement
tic;

%% read the video
myvideo = VideoReader( inputVideoFileName );
if (settings.numFrames < 0)
    settings.numFrames = round(myvideo.Duration*myvideo.FrameRate); %% here you can optionally use less frames, to speedup the processing for testing
end

%% initialize the input video
inputVideo = zeros(round(myvideo.Height*settings.scaleFactor), round(myvideo.Width*settings.scaleFactor), round(settings.numFrames));
disp('Loading video file ...');
for i=1:settings.numFrames
    
    %% read the current frame
    currentFrame = read( myvideo, settings.startFrame+i-1);
    
    %% check if image is rgb and handle accordingly
    if (size(currentFrame,3) == 1)
        inputVideo(:,:,i) = single(imadjust(imresize(im2double(currentFrame), settings.scaleFactor)));
    else
        inputVideo(:,:,i) = single(imadjust(imresize(im2double((rgb2gray( currentFrame ))), settings.scaleFactor)));
    end
    
    %% invert the intensity
    if (settings.invertIntensity == true)
        inputVideo(:,:,i) = 1-inputVideo(:,:,i);
    end
    
    %% display status
    if (mod(i,10) == 0)
        disp(['Loading video file ' num2str(round(100*i/settings.numFrames)) '% ...']);
    end
end

%% calculate the average image of the video
meanImage = single(mean(inputVideo, 3));
imageCenter = size(meanImage) / 2;

%% iterate over all video frames and extract the potential fish locations
results = struct();
%parfor i=1:settings.numFrames %% started at 1
for i=1:settings.numFrames
    
    %% get the current frame and subtract the mean image
    %% the image is smoothed and the binary regions serve as detected objects
    if (settings.usePartialMeanImage == true)
        currentFrame = double(squeeze(inputVideo(:,:,i)) - mean(squeeze(inputVideo(:,:,[max(1,(i-settings.partialMeanImageRadius)),min((i+settings.partialMeanImageRadius),settings.numFrames)])), 3));
    else
        currentFrame = double(squeeze(inputVideo(:,:,i)) - meanImage);
    end
    
   % currentFrame = imgaussfilt(currentFrame, settings.gaussianSigma);
    currentFrame(currentFrame < 0) = 0;
    currentFrame = imfilter(currentFrame, filterMask);
    
    figure(1); imagesc(currentFrame);
    
    if (settings.useLocalAdaptiveThreshold == true)
        currentFrameThresholded = adaptivethreshold(currentFrame,settings.localAdaptiveThresholdWindowSize, 0, 0);
    else
        currentFrameThresholded = im2bw(currentFrame, max(settings.minimumIntensityThreshold, graythresh(currentFrame)));  %% graythresh = 0.03 for hyperactivity screen! %% max(0.05, graythresh(currentFrame)/2))         %%currentFrame >= movementThreshold;
    end
    
    currentFrameWatershed = double(watershed(1-currentFrame));
    if (settings.useWatershedBasedSplitting == true)
        currentFrameThresholded = currentFrameWatershed.*double(currentFrameThresholded) > 0;
    end
    currentFrameRegionProps = regionprops(currentFrameThresholded, currentFrame, 'WeightedCentroid', 'Area', 'MaxIntensity');
    
    %% result matrix is filled with the extracted regionprops
    resultMatrix = [];
    for j=1:length(currentFrameRegionProps)
        
        %% extract maximum intensity of a region surrounding the current centroid to supress fragments that do not contain the head region
        rangeY = max(1,round(currentFrameRegionProps(j).WeightedCentroid(1)-settings.objectRadius)):min(round(currentFrameRegionProps(j).WeightedCentroid(1)+settings.objectRadius), size(currentFrame,2));
        rangeX = max(1,round(currentFrameRegionProps(j).WeightedCentroid(2)-settings.objectRadius)):min(round(currentFrameRegionProps(j).WeightedCentroid(2)+settings.objectRadius), size(currentFrame,1));
        neighborhoodMaxIntensity = max(max(currentFrame(rangeX, rangeY)));
        
        watershedRegionMaxIntensity = max(max(currentFrame(currentFrameWatershed == currentFrameWatershed(round(currentFrameRegionProps(j).WeightedCentroid(2)), round(currentFrameRegionProps(j).WeightedCentroid(1))))));
        intensityRatio = watershedRegionMaxIntensity / neighborhoodMaxIntensity;
        
        if (intensityRatio > 0.5 && currentFrameRegionProps(j).Area >= settings.minArea)  %% currentFrameRegionProps(j).MaxIntensity >= watershedRegionMaxIntensity
            resultMatrix = [resultMatrix; j, currentFrameRegionProps(j).Area, currentFrameRegionProps(j).WeightedCentroid(1), currentFrameRegionProps(j).WeightedCentroid(2), inf, 0];
        end
    end
    
    %% save the results to the data structure
    if (settings.debugFigures == true)
        results(i).DifferenceImage = currentFrame;
        results(i).BinaryImage = currentFrameThresholded;
    end
    if (~isempty(resultMatrix))
        results(i).RegionProps = sortrows(resultMatrix, -2);
    else
        results(i).RegionProps = [];
    end
    
    %% display progress message
    disp(['Finished processing frame ' num2str(i)]);
end

if (settings.debugFigures == true)
    
    meanIntensity = single(squeeze(mean(mean(inputVideo, 2),1)));
    
end
%     figure(2);
%     subplot(1,2,1);
%     plot(1:length(meanIntensity), meanIntensity, '-r');
%     ylabel('Mean Intensity (a.u.)');
%     xlabel('Frame Number');
%
%     subplot(1,2,2);
%     plot(1:length(thresholdValues), thresholdValues, '-r');
%     ylabel('Threshold Values');
%     xlabel('Frame Number');
% end

%% link the segmented regions
figure(1); hold on;
d_orgs = zeros(0, settings.numFrames-1, 6);

for i=1:(settings.numFrames-1)
    
    currentPositions = results(i).RegionProps;
    nextPositions = results(i+1).RegionProps;
    if (i==1)
        d_orgs(1:size(currentPositions,1),i,:) = currentPositions;
    end
    
    if (~isempty(currentPositions) && ~isempty(nextPositions))
        for j=1:size(currentPositions, 1)
            
            %% get the current position and create a matrix containing the same entry count as the next time point
            currentPosition = ones(size(nextPositions(:,3:4)));
            currentPosition(:,1) = currentPosition(:,1) * currentPositions(j, 3);
            currentPosition(:,2) = currentPosition(:,2) * currentPositions(j, 4);
            
            %% calculate the distances to the next frame and identify the nearest neighbor
            distances = sqrt(sum((nextPositions(:,3:4) - currentPosition).^2, 2));
            [minDistance, minIndex] = min(distances);
            areaRatio = max(currentPositions(j,2), nextPositions(minIndex,2)) / min(currentPositions(j,2), nextPositions(minIndex,2));
            secondNearestNeighborRatio = 0;
            if (length(distances) > 1)
                sortedDistances = sortrows(distances);
                secondNearestNeighborIndex = find(distances == sortedDistances(2));
                secondNearestNeighborRatio = sortedDistances(1)/sortedDistances(2);
            end
            
            %% match position already has a partner, i.e., only assign it if this one is closer
            if (nextPositions(minIndex,6) > 0 && minDistance > nextPositions(minIndex,5))
                continue;
            end
            
            %% only add the match in case the areas are rather similar (otherwise the centroid jumps) and if the maximum intensity is sufficiently large
            if (currentPositions(j,2) > settings.minArea && nextPositions(minIndex,2) > settings.minArea && minDistance <= settings.maxVelocity) %% (areaRatio < maxAreaRatio || secondNearestNeighborRatio < 0.7) &&
                
                if (nextPositions(minIndex,6) > 0 && secondNearestNeighborRatio > 0.6)
                    nextPositions(secondNearestNeighborIndex,5) = sortedDistances(2);
                    nextPositions(secondNearestNeighborIndex,6) = nextPositions(minIndex,6);
                end
                nextPositions(minIndex,5) = minDistance;
                nextPositions(minIndex,6) = j;
            else
                disp(['Skipping object with area: ' num2str(currentPositions(j,2)) ', areaRatio: ' num2str(areaRatio) ', snnr: ' num2str(secondNearestNeighborRatio)]);
            end
        end
    end
    
    %% update the next frame positions
    results(i+1).RegionProps = nextPositions;
    
    %% additionally store all tracked results in matrix format
    d_orgs(1:size(nextPositions,1),i+1,:) = nextPositions;
end

%% display progress message
disp('Finished processing track information ...');

%% calculate the average values per frame
d_orgs_new = zeros(1, size(d_orgs,2), 12);
var_bez = char('numActiveEmbryos', ...
    'numIdleEmbryos', ...
    'totalDistPxPerWell', ...
    'meanDistPxEmbryoNormalized', ...
    'meanDistPxActivityNormalized', ...
    'meanSwimSpeedPxPerSecondEmbryoNormalized', ...
    'meanSwimSpeedPxPerSecondActivityNormalized', ...
    'totalDistMMPerWell', ...
    'meanDistMMEmbryoNormalized', ...
    'meanDistMMActivityNormalized', ...
    'meanSwimSpeedMMPerSecondEmbryoNormalized', ...
    'meanSwimSpeedMMPerSecondActivityNormalized');
for i=1:size(d_orgs,2)
    currentDistances = squeeze(d_orgs(:,i,5));
    currentDistances(isinf(currentDistances)) = 0;
    validIndices = (currentDistances > settings.movementThreshold);
    
    numActiveEmbryos(i) = sum(validIndices);
    totalDistPxPerWell(i) = sum(currentDistances);
    meanDistancePxEmbryoNormalized(i) = mean(currentDistances);
    meanDistancePxActivityNormalized(i) = mean(currentDistances(validIndices));
    
    if (isinf(meanDistancePxEmbryoNormalized(i)) || isnan(meanDistancePxEmbryoNormalized(i)))
        meanDistancePxEmbryoNormalized(i) = 0;
    end
    if (isinf(meanDistancePxActivityNormalized(i)) || isnan(meanDistancePxActivityNormalized(i)))
        meanDistancePxActivityNormalized(i) = 0;
    end
    
    meanSwimSpeedPxPerSecondEmbryoNormalized(i) = meanDistancePxEmbryoNormalized(i) * myvideo.FrameRate;
    meanSwimSpeedPxPerSecondActivityNormalized(i) = meanDistancePxActivityNormalized(i) * myvideo.FrameRate;
    
    d_orgs_new(1, i, :) = [numActiveEmbryos(i), ...
        settings.numObjects - numActiveEmbryos(i), ...
        totalDistPxPerWell(i), ...
        meanDistancePxEmbryoNormalized(i), ...
        meanDistancePxActivityNormalized(i), ...
        meanSwimSpeedPxPerSecondEmbryoNormalized(i), ...
        meanSwimSpeedPxPerSecondActivityNormalized(i), ...
        totalDistPxPerWell(i), ...
        meanDistancePxEmbryoNormalized(i), ...
        meanDistancePxActivityNormalized(i), ...
        meanSwimSpeedPxPerSecondEmbryoNormalized(i), ...
        meanSwimSpeedPxPerSecondActivityNormalized(i)];
    
    d_orgs_new(1, i, (end-4):end) = d_orgs_new(1, i, (end-4):end) * settings.pixelToMMFactor;
end

%% the average result values
dorgbez = char('activityPercentage', ...
    'idlePercentage', ...
    'activeEmbryosPercentage', ...
    'idleEmbryosPercentage', ...
    'totalDistancePxPerWell', ...
    'totalDistancePxEmbryoNormalized', ...
    'totalDistancePxActivityNormalized', ...
    'meanSwimSpeedPxEmbryoNormalized', ...
    'meanSwimSpeedPxActivityNormalized', ...
    'totalDistanceMMPerWell', ...
    'totalDistanceMMEmbryoNormalized', ...
    'totalDistanceMMActivityNormalized', ...
    'meanSwimSpeedMMEmbryoNormalized', ...
    'meanSwimSpeedMMActivityNormalized');
d_org_new = zeros(1,14);
d_org_new(1,1) = 100 * sum(squeeze(d_orgs_new(1, :, 1)) > 0) / settings.numFrames;
d_org_new(1,2) = 100 - d_org_new(1,1);
d_org_new(1,3) = 100 * mean(squeeze(d_orgs_new(1, :, 1))) / settings.numObjects;
d_org_new(1,4) = 100 - d_org_new(1,3);
d_org_new(1,5) = sum(squeeze(d_orgs_new(1, :, 3)));
d_org_new(1,6) = sum(squeeze(d_orgs_new(1, :, 4)));
d_org_new(1,7) = sum(squeeze(d_orgs_new(1, :, 5)));
d_org_new(1,8) = mean(squeeze(d_orgs_new(1, :, 6)));
d_org_new(1,9) = mean(squeeze(d_orgs_new(1, :, 7)));
d_org_new(1,10) = sum(squeeze(d_orgs_new(1, :, 8)));
d_org_new(1,11) = sum(squeeze(d_orgs_new(1, :, 9)));
d_org_new(1,12) = sum(squeeze(d_orgs_new(1, :, 10)));
d_org_new(1,13) = mean(squeeze(d_orgs_new(1, :, 11)));
d_org_new(1,14) = mean(squeeze(d_orgs_new(1, :, 12)));
code = 1;

d_orgs_old = d_orgs;
d_orgs = d_orgs_new;
d_org = d_org_new;
save(strrep(settings.resultAviFileName, '.mp4', '.prjz'), '-mat', 'd_orgs', 'd_org', 'var_bez', 'dorgbez', 'code');
d_orgs = d_orgs_old;

%% averaged result values
d_orgs(isinf(d_orgs(:))) = 0;
d_orgs(:,:,5) = d_orgs(:,:,5) / settings.scaleFactor;
averageResults = [];
averageResults.activePercentage = 100 * sum(sum((d_orgs(:,:,5) > settings.movementThreshold))) / settings.numObjects / settings.numFrames; %% returns the percentage of frames the embryos were active
averageResults.idlePercentage = 100 - averageResults.activePercentage; %% returns the percentage of frames the embryos were active
distances = squeeze(d_orgs(:,:,5));
distances(isinf(distances)) = 0;
averageResults.totalDistancePixel = sum(sum((distances))) / settings.numObjects; % returns the average distance traveled in pixels which are scaled by the scaleFactor
averageResults.totalDistanceMM = averageResults.totalDistancePixel * settings.pixelToMMFactor; % returns the average distance traveled in mm
averageResults.swimSpeedPixelPerSecond = mean(distances(distances > settings.movementThreshold)) * myvideo.FrameRate;    %% mean instantaneous velocity per active swimming (pixel/s)
if (isnan(averageResults.swimSpeedPixelPerSecond))
    averageResults.swimSpeedPixelPerSecond = 0;
end
averageResults.swimSpeedMMPerSecond = averageResults.swimSpeedPixelPerSecond * settings.pixelToMMFactor;  %% mean instantaneous velocity per active swimming (mm/s)

disp('Finished calculating average statistics ...');

%% show starting message
disp('Starting to create results plot ...');

%% open new figure for plotting
fh = figure;
set(fh, 'units','pixel', 'OuterPosition', [0    0    1024    512]);
colordef black;
set(fh,'Color', [0,0,0]);

%% open subplot for the mean image
subplot(1,2,1);
imagesc(meanImage); hold on;
axis equal;
subplot(1,2,2);
hold on;

%% plot the movement lines for all embryos
colormap gray;
for i=2:settings.numFrames
    subplot(1,2,1);
    for j=1:size(d_orgs,1)
        
        if (d_orgs(j,i,3) ~= 0 && d_orgs(j,i,4) ~=0)
            plot(d_orgs(j,i,3), d_orgs(j,i,4), '.r');
        end
        
        if (round(d_orgs(j,i,6)) > 0 && round(d_orgs(j,i,6)) <= size(d_orgs,1))
            predecessor = round(d_orgs(j,i,6));
            if (d_orgs(predecessor, i-1,3) == 0 || d_orgs(predecessor, i-1,4) == 0 || d_orgs(j, i,3) == 0 || d_orgs(j, i,4) == 0)
                continue;
            end
            
            plot([d_orgs(predecessor, i-1,3), d_orgs(j, i,3)], [d_orgs(predecessor, i-1,4), d_orgs(j, i,4)], '-.m');
        end
    end
    
    subplot(1,2,2);
    plot([i-1,i], [sum(d_orgs(:,i-1,5))/settings.numObjects, sum(d_orgs(:,i,5))/settings.numObjects], '-c', 'LineWidth', 2);
end

%% set axis labels for the average displacement plot
subplot(1,2,2);
axis([-1, settings.numFrames, 0, 20]);
xlabel('Frame number');
ylabel('Average displacement (px)');
legend('Average displacement (px)');

%% get the current frame
F = getframe(fh);
[resultImage, ~] = frame2im(F);
imwrite(resultImage, strrep(settings.resultAviFileName, '.mp4', '.png'));

%% close the figure
%close(fh);
disp('Finished creating result plot ...');

%% result movie generation (optional)
MotilityAnalysisDebugVideo;
end