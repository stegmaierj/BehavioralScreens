function larvaModels = IdentifyLarvaModels(labelImage, rawImage, maskImage, numSegments, segmentRadius, roiRadius, larvaLength)

%% enable/disable debug flag
debugFigures = false;

%% set default parameters
if nargin < 4
    numSegments = 2;
end
if nargin < 5
    segmentRadius = 10;
end
if nargin < 6
    roiRadius = 10;
end
if nargin < 7
    larvaLength = 35;
end

%% maximum angle of the line segments
maximumAngle = pi/2;

%% identify larva segments
currentRegionProps = regionprops(labelImage>0, rawImage, 'Centroid', 'Area', 'PixelIdxList', 'MeanIntensity', 'Eccentricity', 'Solidity', 'ConvexArea', 'EquivDiameter');
endPointsImage = bwmorph(labelImage > 0, 'endpoints') .* labelImage;
currentLarva = 1;
larvaModels = struct();

%% calculate the difference of gaussian image for more accurate head and first segment localization
dogImage = imgaussfilt(rawImage, 5) - imgaussfilt(rawImage, 3);
dogImage = imgaussfilt((1 - (dogImage .* maskImage)), 2);

%% iterate over all regions and identify the head locations
for r=1:length(currentRegionProps)
    
    %% identify the number of potential larva in the roi
    numLarva = max(1, round(currentRegionProps(r).Area / larvaLength));
    
    %% identify the potential end points
    [currentEndPointsY, currentEndPointsX] = ind2sub(size(endPointsImage), find(endPointsImage == r));
    
    %% check if at least one end point exists
    if (length(currentEndPointsX) > 1)

        %% perform clustering to identify the potential head regions
        currentLinkage = linkage([currentEndPointsX, currentEndPointsY], 'ward', 'euclidean');
        currentClustering = cluster(currentLinkage, 'Cutoff',7,'Criterion','distance');
        clusterIndices = unique(currentClustering);
        clusterMeans = [];
        clusterIntensities = [];
        for i=clusterIndices'
            clusterMeans = [clusterMeans; round(mean([currentEndPointsX(currentClustering == i), currentEndPointsY(currentClustering == i)],1))];
            clusterIntensities = [clusterIntensities; rawImage(clusterMeans(end,2), clusterMeans(end,1))];
        end
        
        %% sort the clustered intensities
        [sortedIntensities, sortedIndices] = sort(clusterIntensities, 'descend');
    else
        clusterMeans = [currentEndPointsX, currentEndPointsY];
        sortedIndices = 1;
    end
    
    %% extract the larvae
    for c=1:min(numLarva, length(sortedIndices))
        
        %% initialize the larva with the identified head position
        endPointIndex = sortedIndices(c);
        initialHeadLocation = clusterMeans(endPointIndex,:);
        
        %% perform more accurate localization of the head position using the smoothed dog image
        headSearchRadius = segmentRadius/2;
        rangeX = max(1, initialHeadLocation(1)-headSearchRadius):min(size(rawImage,1), initialHeadLocation(1)+headSearchRadius);
        rangeY = max(1, initialHeadLocation(2)-headSearchRadius):min(size(rawImage,2), initialHeadLocation(2)+headSearchRadius);
        headSearchRegion = dogImage(rangeY, rangeX);
        [maxValues, maxIndices] = max(headSearchRegion(:));
        [xpos, ypos] = ind2sub(size(headSearchRegion), maxIndices);
        correctedHeadLocation = [initialHeadLocation(1) + (ypos - ((length(rangeY)-1)/2) - 1), initialHeadLocation(2) + (xpos - ((length(rangeX)-1)/2) - 1)];
        larvaModel = [correctedHeadLocation];
        directionAngles = [];
        currentPoint = correctedHeadLocation; %clusterMeans(endPointIndex,:);

        %% for debug visualization
        if (debugFigures == true)
            figure(2);
            imagesc(rawImage); hold on;
            plot(currentPoint(1), currentPoint(2), '*r');
        end
        
        %% try to extract the respective segments starting from head and moving towards the tail
        for i=1:numSegments
            
            %% use a larger step at the first point to estimate the initial direction correctly
            if (i==1)
                currentSegmentRadius = segmentRadius * 1;
            else
                currentSegmentRadius = segmentRadius;
            end
            
            %% specify the search range based on the segment radius
            rangeX = max(1, currentPoint(1)-currentSegmentRadius):min(size(rawImage,1), currentPoint(1)+currentSegmentRadius);
            rangeY = max(1, currentPoint(2)-currentSegmentRadius):min(size(rawImage,2), currentPoint(2)+currentSegmentRadius);
            
            %% get regional intensities
            minRegionIntensity = min(min(rawImage(rangeY, rangeX)));
            maxRegionIntensity = max(max(rawImage(rangeY, rangeX)));
            
            %% initialize the potential next point and extremum search variables
            nextPoint = currentPoint;
            maxValue = 0;
            maxAngle = 0;
            angleOverView = zeros(size(rawImage));
            decisionFSMDImage = zeros(size(rawImage));
            
            %% iterate over the region and identify the most likely move
            for j=rangeX
                for k=rangeY
                    potentialPoint = [j,k];
                    if ((norm(potentialPoint-currentPoint) <= currentSegmentRadius+1 && ...
                        norm(potentialPoint-currentPoint) > (currentSegmentRadius-1)) || ...
                        i == 1)
                        
                        %% calculate average intensity along the line
                        [intX, intY] = bresenham(larvaModel(i,1), larvaModel(i,2), j, k);
                        averageIntensitySkeleton = 0;
                        averageIntensityRawImage = 0;
                        for l=1:length(intX)
                            averageIntensitySkeleton = averageIntensitySkeleton + (dogImage(intY(l), intX(l))); %% (labelImage(intY(l), intX(l)) > 0)
                            averageIntensityRawImage = averageIntensityRawImage + rawImage(intY(l), intX(l));   %% rawImage(intY(l), intX(l))
                        end
                        %averageIntensitySkeleton = averageIntensitySkeleton / length(intX);
                        %averageIntensityRawImage = averageIntensityRawImage / length(intX);
                        %averageIntensityRawImage = rawImage(k,j);
                        
                        % calculate the FSMD values
                        if (i > 1)
                            previousPointDirection = (larvaModel(i,:)-larvaModel(i-1,:)) / norm(larvaModel(i,:)-larvaModel(i-1,:));
                            nextPointDirection = (potentialPoint-larvaModel(i,:))/norm(potentialPoint-larvaModel(i,:));
                            angle = real(acos(dot(previousPointDirection, nextPointDirection)));
                            
                            if (debugFigures == true)
                                angleOverView(k,j) = angle;
                            end
                            
                            if (isnan(angle))%(angle > maximumAngle || isnan(angle))
                                continue;
                            end
                            
                            %decisionFSMD = (trapmf(angle, [0,0,pi/2,pi])) * max(averageIntensitySkeleton, trapmf(averageIntensityRawImage, [minRegionIntensity, maxRegionIntensity, maxRegionIntensity, maxRegionIntensity]));
                            decisionFSMD = (trapmf(angle, [0,0,pi/2,pi])) * averageIntensityRawImage;
                        else
                            decisionFSMD = averageIntensitySkeleton; %dogImage(k, j); %labelImage(k,j) > 0; %averageIntensitySkeleton;%labelImage(k,j);
                        end
                        
                        if (~isreal(decisionFSMD))
                            test = 1;
                        end
                        
                        decisionFSMDImage(k,j) = decisionFSMD;
                        
                        if (decisionFSMD > maxValue)
                            nextPoint = [j,k];
                            maxValue = decisionFSMD;
                            if (i>1)
                                maxAngle = angle;
                            end
                        end
                    end
                end
            end
            
            %% for debug visualizations
            if (debugFigures == 1)
                figure(3);
                imagesc(rawImage); hold on;
                plot(initialHeadLocation(1), initialHeadLocation(2), '.c');
                plot([currentPoint(1), nextPoint(1)], [currentPoint(2), nextPoint(2)], '-or');
            end
            
            %% only add the segment if it's at least larger than zero
            if (maxValue > 0)
                larvaModel = [larvaModel; nextPoint];
                if (i>1)
                    directionAngles = [directionAngles; maxAngle];
                end
                currentPoint = nextPoint;
            else
                break;
            end
        end
        
        %% check if the larva exits the region of interest and discard it in this case
        validLarva = true;
        larvaDescriptor = zeros(size(larvaModel,1), 1);
        for i=1:size(larvaModel,1)
            larvaDescriptor(i) = rawImage(larvaModel(i,1), larvaModel(i,2));
            if (maskImage(larvaModel(i,1), larvaModel(i,2)) == 0)
                validLarva = false;
                break;
            end
        end
        
        %% add the larva to the results data structure
        if (validLarva == true)
            %% copy the larva model
            larvaModels(currentLarva).larvaModel = larvaModel;
            
            %% calculate the direction angles
            directionAngles = [];
            relativeAngles = [];
            currentLength = 0;
            for i=2:size(larvaModel,1)
                currentDirection = (larvaModel(i,:)-larvaModel(i-1,:));
                if i<size(larvaModel,1)
                    nextDirection = (larvaModel(i+1,:)-larvaModel(i,:));
                    relativeAngles = [relativeAngles; radtodeg(real(acos(dot((nextDirection/norm(nextDirection)), (currentDirection/norm(currentDirection))))))];
                end
                
                directionLength = norm(larvaModel(i,:)-larvaModel(i-1,:));
                currentDirection = currentDirection / directionLength;
                currentLength = currentLength + directionLength;
                directionAngles = [directionAngles; radtodeg(acos(dot(currentDirection, [1,0])))];
            end
            angleChange = directionAngles - circshift(directionAngles,-1);
            
            %% set features of the current larva
            larvaModels(currentLarva).directionAngles = directionAngles;
            larvaModels(currentLarva).relativeAngles = relativeAngles;
            larvaModels(currentLarva).relativeAngleSum = sum(relativeAngles);
            larvaModels(currentLarva).angleChange = angleChange(1:end-1);
            larvaModels(currentLarva).angleChangeSum = sum(angleChange(1:end-1));
            larvaModels(currentLarva).length = currentLength;
            larvaModels(currentLarva).area = currentRegionProps(r).Area;
            larvaModels(currentLarva).eccentricity = currentRegionProps(r).Eccentricity;
            larvaModels(currentLarva).solidity = currentRegionProps(r).Solidity;
            larvaModels(currentLarva).convexArea = currentRegionProps(r).ConvexArea;
            larvaModels(currentLarva).equivDiameter = currentRegionProps(r).EquivDiameter;
            currentLarva = currentLarva + 1;
        end
        if (debugFigures == true)
            plot(larvaModel(:,1), larvaModel(:,2), '-r');
            test = 1;
        end
    end
end
end