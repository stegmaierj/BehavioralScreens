function circleROIImage = ExtractWellMask(meanImage, roiRadiusPadding)

    openedMeanImage = imopen(im2bw(meanImage, quantile(meanImage(:), 0.3)), strel('disk', 3));
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
                circleROIImage(j,i) = 1;
            end
        end
    end

end