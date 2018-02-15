
function [] = TouchResponse(inputDir, outputDir, markMode)
%% specify the input and output folder
inputFiles = dir(inputDir);

%% if true, only the first frames are loaded for landmark placement.
%markMode = true;
debugFigures = true;
imagePadding = 32;
templateRadiusEmbryoMax = 5;
templateRadiusMax = 15;
templateRadiusNeedleMax = 5;
resultVideoFrameStride = 20;
startFrame = 1;

%% process all files in the folder
for f=3:length(inputFiles)
    
    %% continue if it's not a valid video file
    if (isempty(strfind(inputFiles(f).name, '.avi')))
        continue;
    end
    
    %% specify the file names for the current input and output files    
    currentFile = [inputDir inputFiles(f).name];
    currentResultFile = strrep([outputDir inputFiles(f).name], '.avi', '.mat');
    currentResultVideoFile = strrep([outputDir inputFiles(f).name], '.avi', '_Result.avi');
    
    %% open the current video
    v = VideoReader(currentFile);
    numFrames = floor(v.FrameRate * v.Duration)-10;
    
    %% initialize the needle and embryo coordinates
    coordsNeedle = zeros(numFrames, 2);
    coordsEmbryo = zeros(numFrames, 2);
    areaNeedle = zeros(numFrames, 1);
    areaEmbryo = zeros(numFrames, 1);
    
    %% open the result video if debug videos should be generated
    if (debugFigures == true)
        videoWriter = VideoWriter(currentResultVideoFile);
        open(videoWriter);
    end
    
    %% specify the start and end frames
    if (markMode == true)
        endFrame = startFrame;
    else
        endFrame = numFrames;
    end
    
    %% loop through all frames and track the needle/embryo
    currentFrame = 1;
    for i=startFrame:endFrame
        
        %% load the current video frame
        currentImage = padarray(imgaussfilt(double(rgb2gray(read(v, i))) / 255, 1), [imagePadding, imagePadding], 'replicate');
        %currentBwImage = zeros(size(currentImage));
        currentBwImage = 1-im2bw(currentImage, graythresh(currentImage));
       
        %% if first frame, open click positions
        if (startFrame == i && (markMode == true || ~exist(currentResultFile, 'file')))
            
            %% open the figure and get coordinates of the needle
            figure(1); clf;
            imagesc(currentImage); title('Mark center of the needle...');
            coordsNeedle(currentFrame,:) = round(ginput(1));
            
            figure(1); title('Mark center of the embryo head...');
            coordsEmbryo(currentFrame,:) = round(ginput(1));
            save(currentResultFile, '-mat', 'coordsNeedle', 'coordsEmbryo');
            break;
        elseif (startFrame == i && exist(currentResultFile, 'file'))
            load(currentResultFile);
            touchTimePoint = 0;
        end
        
        if (startFrame == i)
            initialImage = currentImage;
            lastImage = currentImage;
        end
        
        %% index of the last frame
        lastFrame = max(currentFrame-1, 1);
        
        distanceNeedleEmbryo = norm(coordsNeedle(lastFrame, :) - coordsEmbryo(lastFrame, :));
        templateRadiusNeedle = min(round(distanceNeedleEmbryo/2), templateRadiusMax);
        
        %% Detect and track the Needle
        rangeX = max(1, coordsNeedle(lastFrame, 1)-templateRadiusNeedleMax):min(size(currentImage,1), coordsNeedle(lastFrame, 1)+templateRadiusNeedleMax);
        rangeY = max(1, coordsNeedle(lastFrame, 2)-templateRadiusNeedleMax):min(size(currentImage,2), coordsNeedle(lastFrame, 2)+templateRadiusNeedleMax);
        validPositionsSearchRange = [];
        for j=rangeX
            for k=rangeY
                validPositionsSearchRange = [validPositionsSearchRange; j, k];
                
%                 if (norm([j,k] - coordsNeedle(lastFrame, :)) < templateRadiusNeedleMax)
%                     currentBwImage(k,j) = 1;
%                 end
            end
        end
        
        %% extract descriptor for the needle in the first frame and search for it in the current frame using the last known position
        [featuresNeedle,~] = extractFeatures(initialImage, (coordsNeedle(1, :)), 'Method', 'block');
        [neighborsNeedle,~] = extractFeatures(currentImage, (validPositionsSearchRange), 'Method', 'block');
        [indexPairs, matchMetric] = matchFeatures(featuresNeedle, neighborsNeedle, 'Unique', true, 'MatchThreshold', 100, 'MaxRatio', 1.0);
        
        %% set the updated needle position
        if (~isempty(indexPairs))
            coordsNeedle(currentFrame, :) = validPositionsSearchRange(indexPairs(2),:);
        else
            coordsNeedle(currentFrame, :) = coordsNeedle(currentFrame-1, :);
        end
        
        currentImageTemp = currentImage;
        currentImage(rangeY, rangeX) = 1;
        
        %% Detect and track the Embryo
        rangeX = max(1, coordsEmbryo(lastFrame, 1)-templateRadiusMax):min(size(currentImage,1), coordsEmbryo(lastFrame, 1)+templateRadiusMax);
        rangeY = max(1, coordsEmbryo(lastFrame, 2)-templateRadiusMax):min(size(currentImage,2), coordsEmbryo(lastFrame, 2)+templateRadiusMax);
        
        %% initialize the embryo mask
%         embryoMask = zeros(size(currentBwImage));
%         embryoMask(rangeY, rangeX) = 1;
%         embryoThreshold = graythresh(currentImage(embryoMask > 0));
%         embryoMask = (embryoMask > 0) .* (1-im2bw(currentImage .* (embryoMask > 0), embryoThreshold));
%         currentBwImage = embryoMask | currentBwImage;
%         embryoMask = imopen(embryoMask, strel('disk', 5));

        %% identify the potential linking positions
        validPositionsEmbryo = [];
        validPositionsSearchRange = [];
        backupPositions = [];
        maxPosition = [0,0];
        maxValue = 0;
        for j=rangeX
            for k=rangeY
                
                %% the search range in the current frame
                if (currentBwImage(k, j) > 0 && norm([j,k] - coordsEmbryo(lastFrame, :)) < templateRadiusMax)
                    validPositionsSearchRange = [validPositionsSearchRange; j, k];
                end
                
                %% if no positions satisfy the conditions above, save some backup positions
                backupPositions = [backupPositions; j, k];
                
                %% positions located on the embryo head in the last frame
                if (norm([j,k] - coordsEmbryo(lastFrame, :)) < templateRadiusEmbryoMax)
                    validPositionsEmbryo = [validPositionsEmbryo; j, k];
                end
            end
        end
        
        %% switch to backup positions if no valid locations are found
        if (isempty(validPositionsSearchRange))
            validPositionsSearchRange = backupPositions;
        end
                
        %% extract descriptor of the embryo head and search for it in the current frame
        [featuresEmbryo1, validPointsEmbryo1] = extractFeatures(lastImage, SURFPoints(validPositionsEmbryo), 'Method', 'SURF', 'SURFSize', 128); %, 'BlockSize', blockSize, 'SURFSize', 128);
        if (startFrame == i)
            [featuresEmbryo2, validPointsEmbryo2] = extractFeatures(initialImage, SURFPoints(validPositionsEmbryo), 'Method', 'SURF', 'SURFSize', 128);%, 'SURFSize', 128);%, 'BlockSize', blockSize, 'SURFSize', 128);
            initialPositionsEmbryo = validPositionsEmbryo;
        end
        [neighborsEmbryo1, validPointsNeighborsEmbryo1] = extractFeatures(currentImage, SURFPoints(validPositionsSearchRange), 'Method', 'SURF', 'SURFSize', 128);%, 'BlockSize', blockSize, 'SURFSize', 128);
        %[neighborsEmbryo2, validPointsNeighborsEmbryo2] = extractFeatures(currentImage, SURFPoints(validPositionsSearchRange), 'Method', 'SURF', 'SURFSize', 128);%, 'SURFSize', 128);%, 'BlockSize', blockSize, 'SURFSize', 128);
        [indexPairs1, matchMetric1] = matchFeatures(featuresEmbryo1, neighborsEmbryo1, 'MatchThreshold', 20.0, 'MaxRatio', 1.0, 'Unique', true);
        [indexPairs2, matchMetric2] = matchFeatures(featuresEmbryo2, neighborsEmbryo1, 'MatchThreshold', 10.0, 'MaxRatio', 0.8, 'Unique', true);
                
        %% set the new embryo position and use old position in case no match was found
        normalizationFactor = 0;
        lastPosition = coordsEmbryo(lastFrame, :);
        coordsEmbryo(currentFrame, :) = 0;
        if (~isempty(indexPairs2) || ~isempty(indexPairs1))
            for j=1:size(indexPairs2,1)
                currentDistance = norm(initialPositionsEmbryo(indexPairs2(j, 1),:) - coordsEmbryo(1, :));
                currentWeight = normpdf(currentDistance, 0, 1);
                normalizationFactor = normalizationFactor + currentWeight;
                coordsEmbryo(currentFrame, :) = coordsEmbryo(currentFrame, :) + currentWeight * validPositionsSearchRange(indexPairs2(j, 2),:);
            end
        %	coordsEmbryo(currentFrame, :) = round(coordsEmbryo(currentFrame, :) / normalizationFactor);
        %elseif (~isempty(indexPairs1))
            %disp('Using initial location for matching');

            for j=1:size(indexPairs1,1)
                currentDistance = norm(validPositionsEmbryo(indexPairs1(j, 1),:) - lastPosition);
                currentWeight = normpdf(currentDistance, 0, 1);
                normalizationFactor = normalizationFactor + currentWeight;
                coordsEmbryo(currentFrame, :) = coordsEmbryo(currentFrame, :) + currentWeight * validPositionsSearchRange(indexPairs1(j, 2),:);
            end
            coordsEmbryo(currentFrame, :) = round(coordsEmbryo(currentFrame, :) / normalizationFactor);
            %disp('Using current location for matching');
        else
            coordsEmbryo(currentFrame, :) = lastPosition;
            %disp('No match found, not moving...');
        end
% %         
%         if (~isempty(indexPairs1) || ~isempty(indexPairs2))
%             figure(2);
%             showMatchedFeatures(lastImage, currentImage, [validPositionsEmbryo(indexPairs1(:, 1),:); validPositionsEmbryo(indexPairs2(:, 1),:)], [validPositionsSearchRange(indexPairs1(:, 2),:); validPositionsSearchRange(indexPairs2(:, 2),:)]);
%         end
        
        %% restore the current image such that the needle roi is not blanked
        currentImage = currentImageTemp;
        
        %% create a label image containing only the embryo and the needle
        currentLabelImage = bwlabel(currentBwImage > 0);
        currentRegionProps = regionprops(currentLabelImage, 'Area');

        %% identify the labels of the embryo and the needle
        needleLabel = currentLabelImage(coordsNeedle(currentFrame, 2), coordsNeedle(currentFrame, 1));
        embryoLabel = currentLabelImage(coordsEmbryo(currentFrame, 2), coordsEmbryo(currentFrame, 1));
        
        %% save region areas in case they are not touching yet
        if (touchTimePoint == 0 && needleLabel > 0)
            areaNeedle(currentFrame, :) = currentRegionProps(needleLabel).Area;
        end
        if (touchTimePoint == 0 && embryoLabel > 0)
            areaEmbryo(currentFrame, :) = currentRegionProps(embryoLabel).Area;
        end
        
        %% set the touch time point to the first frame the two segments of embryo and needle are intersecting
        if (touchTimePoint == 0 && (needleLabel == embryoLabel) && needleLabel > 0 && embryoLabel > 0)
            touchTimePoint = currentFrame;
        end
                
        %% calculate embryo and needle speed
        distanceNeedleEmbryo = sqrt(sum((coordsNeedle(1:currentFrame,:) - coordsEmbryo(1:currentFrame,:)).^2,2));
        speedNeedle = sqrt(sum((diff(coordsNeedle , 1, 1).^2), 2));
        speedEmbryo = sqrt(sum((diff(coordsEmbryo, 1, 1).^2), 2));
        speedNeedle(currentFrame:end) = 0;
        speedEmbryo(currentFrame:end) = 0;
        
        %% generate debug video
        if (debugFigures == true && mod(currentFrame, resultVideoFrameStride) == 0)
            fh = figure(1); set(gcf, 'OuterPosition', [-1626         423        1177         567]);
            clf;
            subplot(3,2,[1,3,5]);
            colormap gray;
            imagesc(currentImage); hold on;
%             currentLabelImage = bwlabel(embryoMask);
%             embryoLabel = currentLabelImage(coordsEmbryo(currentFrame, 2), coordsEmbryo(currentFrame, 1));
%             currentBoundary = bwboundaries(currentLabelImage == embryoLabel);
%             for j=1:length(currentBoundary)
%                 plot(currentBoundary{j}(:,2), currentBoundary{j}(:,1), 'g', 'LineWidth', 1)
%             end
            
            plot(coordsEmbryo(currentFrame,1), coordsEmbryo(currentFrame,2), 'og');
            plot(coordsNeedle(currentFrame,1), coordsNeedle(currentFrame,2), 'ob');
            plot(coordsNeedle(1:currentFrame,1), coordsNeedle(1:currentFrame,2), '-b');
            plot(coordsEmbryo(1:currentFrame,1), coordsEmbryo(1:currentFrame,2), '-g');
            axis equal;
            axis tight
            
            subplot(3,2,2); hold on;
            plot(distanceNeedleEmbryo, '-b');
            if (touchTimePoint > 0)
                plot([touchTimePoint, touchTimePoint], [0, max(distanceNeedleEmbryo)+10], '-k');
                legend('Distance Embryo-Needle', 'Touch Time Point', 'Location', 'NorthWest');
            else
                legend('Distance Embryo-Needle', 'Location', 'NorthWest');
            end
            xlabel('Frame Number');
            ylabel('Distance');
            set(gca, 'XLim', [0, numFrames]);
            set(gca, 'YLim', [0, max(distanceNeedleEmbryo)+10]);
            
            subplot(3,2,4); hold on;
            plot(areaNeedle(1:currentFrame), '-b');
            plot(areaEmbryo(1:currentFrame), '-g');
            if (touchTimePoint > 0)
                plot([touchTimePoint, touchTimePoint], [0, max(max(areaEmbryo, areaNeedle))+10], '-k');
                legend('Needle', 'Embryo', 'Touch Time Point', 'Location', 'NorthWest');
            else
                legend('Needle', 'Embryo', 'Location', 'NorthWest');
            end
            xlabel('Frame Number');
            ylabel('Area');
            set(gca, 'XLim', [0, numFrames]);
            set(gca, 'YLim', [0, max(max(areaEmbryo, areaNeedle))+10]);
            
            
            subplot(3,2,6); hold on;
            plot((speedNeedle(1:min(length(speedNeedle), currentFrame))), '-b');
            plot((speedEmbryo(1:min(length(speedNeedle), currentFrame))), '-g');
            if (touchTimePoint > 0)
                plot([touchTimePoint, touchTimePoint], [0, max(max(speedEmbryo, speedNeedle))+1], '-k');
                legend('Needle', 'Embryo', 'Touch Time Point', 'Location', 'NorthWest');
            else
                legend('Needle', 'Embryo', 'Location', 'NorthWest');
            end
            xlabel('Frame Number');
            ylabel('Speed (Px/Frame)');
            set(gca, 'XLim', [0, numFrames]);
            set(gca, 'YLim', [0, max(max(speedEmbryo, speedNeedle))+1]);
            drawnow;
            
            %% write current frame
            writeVideo(videoWriter, getframe(fh));
        end
        
        %% increase the current frame
        lastImage = currentImage;
        currentFrame = currentFrame+1;
        
        %% display progress
        if (mod(currentFrame,20) == 0)
            disp([num2str(min(100, (100*(currentFrame/numFrames)))) '%']);
        end
    end
    
    %% continue if mark mode is active
    if (markMode == true)
        continue;
    end
    
    %% close the video writer handle
    if (debugFigures == true)
        close(videoWriter);
    end
    
    distanceEmbryo = sum(speedEmbryo);
    distanceNeedle = sum(speedNeedle);
    
    latency = -1;
    if (touchTimePoint > 0)
        speedEmbryoStd = 3*std(speedEmbryo(1:touchTimePoint));
        latency = min(find(speedEmbryo(touchTimePoint:end) > speedEmbryoStd));
    end
    
    %% save the results
    save(currentResultFile, '-mat', 'coordsNeedle', 'coordsEmbryo', 'areaNeedle', 'areaEmbryo', 'speedNeedle', 'speedEmbryo', 'currentFile', 'touchTimePoint', 'latency', 'distanceEmbryo', 'distanceNeedle');
end
end