rawImagesFolder = 'E:\Projects\2016\CBendAnalysis\Data\SingleImages\Video';
inputFolder = 'E:\Projects\2016\CBendAnalysis\Results\item_0018_BinaryThinningImageFilter';
meanImage = imadjust(im2double(imread('E:\Projects\2016\CBendAnalysis\Data\SingleImages\MeanImage.tif')));
rawInputFiles = dir(rawImagesFolder);
inputFiles = dir(inputFolder);

roiRadiusPadding = 4;
minArea = 10;

openedMeanImage = imopen(im2bw(meanImage, graythresh(meanImage)), strel('disk', 3));
meanImageRegionProps = regionprops(openedMeanImage, 'Area', 'Centroid', 'EquivDiameter');

maxIndex = 1;
maxArea = 0;
for i=1:length(meanImageRegionProps)
   if (meanImageRegionProps(i).Area > maxArea)
       maxArea = meanImageRegionProps(i).Area;
       maxIndex = i;
   end
end

circleROIImage = zeros(size(meanImage));
maxRadius = meanImageRegionProps(maxIndex).EquivDiameter / 2 - roiRadiusPadding;
maxCentroid = meanImageRegionProps(maxIndex).Centroid;
for i=1:size(meanImage,1)
    for j=1:size(meanImage,2)
        currentDistance = sqrt(sum((maxCentroid - [i, j]).^2));
        if (currentDistance < maxRadius)
            circleROIImage(i,j) = 1;
        end
    end
end


writerObj = VideoWriter('ResultVideo3.avi', 'MPEG-4');
writerObj.FrameRate = 10;
open(writerObj);

%
 figure; imagesc(circleROIImage);% imagesc(cat(3,circleROIImage, meanImage, meanImage));

fh = figure(1);
lastEccentricity = -1;
for i=3:length(inputFiles)
    
    currentRawImage = im2double(imread([rawImagesFolder filesep rawInputFiles(i).name]));
    currentImage = im2double(imread([inputFolder filesep inputFiles(i).name])) .* circleROIImage;
    
    smoothedRawImage = imgaussfilt(currentRawImage, 5);
        
    tic;
    skeletonImage = bwmorph(currentImage > 0, 'skel');
    toc;
    
    currentLabelImage = bwlabel(skeletonImage);
    currentRegionProps = regionprops(skeletonImage, 'Area', 'PixelIdxList');
    for j=1:length(currentRegionProps)
       if (currentRegionProps(j).Area < minArea)
          currentLabelImage(currentRegionProps(j).PixelIdxList) = 0;
       end
    end    
    
    endPoints = bwmorph(currentLabelImage > 0, 'endpoints');
    subplot(1,2,1); cla;
    set(gca, 'Units', 'normalize', 'Position', [0.025,0.05,0.45,0.925]); 
    imagesc(cat(3,currentLabelImage, endPoints, currentRawImage));
    %imagesc(currentRawImage);
    
    endPointRegionProps = regionprops(endPoints, smoothedRawImage, 'Centroid', 'MeanIntensity');
    
    
    currentRegionProps = regionprops(currentLabelImage > 0, 'Eccentricity');
    averageEccentricity = 0;
    for j=1:length(currentRegionProps)
        averageEccentricity = averageEccentricity + currentRegionProps(j).Eccentricity;
    end
    averageEccentricity = averageEccentricity / length(currentRegionProps);
    if (lastEccentricity < 0)
        lastEccentricity = averageEccentricity;
    end    
    
    subplot(1,2,2);
    if (i==3)
        set(gca, 'Units', 'normalize', 'Position', [0.525,0.05,0.45,0.925]); 
    end
    hold on;
    plot([i-2, i-1], [lastEccentricity, averageEccentricity], '-m');
    ylabel('Eccentricity');
    xlabel('Frame Number');
    axis([0,450,0,1.01]);
    lastEccentricity = averageEccentricity;
    
        writeVideo(writerObj,getframe(fh));
    
end
close(writerObj);