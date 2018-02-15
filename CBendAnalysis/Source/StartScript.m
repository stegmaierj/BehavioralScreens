function [ output_args ] = StartScript( inputFile, generateResultVideo )
%% setup required parameters
settings = struct();
settings.frameRate = 1000; %% the number of frames per second (used to scale latency values to ms)
settings.startFrame = 0; %% skips all frames before that frame to reduce the processing time
settings.numFramesAfterStimulus = 1000; %% the number of frames to consider after the pulse
settings.minPeakHeight = 20; %% minimium angle change in degree. Only c-bends with head angle changes larger than this value are considered a valid c-bend
settings.minPeakProminence = 20; %% minimum relative height of a peak, can be used to supress small peaks between c1 and c2 bends
settings.minTrackletLengthAfterStimulus = 50; %% tracklets that are shorter than this minimum length are discarded, as they potentially do not contain both c-bends
settings.framePadding = 70; %% this amount of frames is used before and after the led-onset, i.e. it should be large enough to cover the desired c-bend behaviors
settings.latencyThreshold = 16; %% head orientation changes above this value are considered as cbend start (16° as described in Burgess07)
settings.roiRadiusPadding = 0; %% used to reduce the circular roi, i.e., the roi mask radius is shrinked by this value in pixels
settings.generateResultVideo = generateResultVideo;

%% setup the input and output paths
[outputDir,outputFile,~] = fileparts(inputFile);

%% initialize dispstat - to be able to overwrite command window messages
dispstat('','init');

%% perform the detection
PerformCBendAnalysis(inputFile, outputDir, settings);

%% perform result computations
PerformResultAnalysis([outputDir '\' outputFile '.mat']);

end