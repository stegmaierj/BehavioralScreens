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

function circleROIImage = ExtractWellMask(meanImage, roiRadiusPadding)

    %% perform opening to remove small structures and disconnect larger connected structures
    openedMeanImage = imopen(im2bw(meanImage, quantile(meanImage(:), 0.3)), strel('disk', 3));
    
    %% identify the connected components
    meanImageRegionProps = regionprops(openedMeanImage, 'Area', 'Centroid', 'EquivDiameter');

    %% find the largest connected component
    maxIndex = 1;
    maxArea = 0;
    for i=1:length(meanImageRegionProps)
       if (meanImageRegionProps(i).Area > maxArea)
           maxArea = meanImageRegionProps(i).Area;
           maxIndex = i;
       end
    end

    %% create the result image containing only the circle of interest
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