
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