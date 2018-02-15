function [ output_args ] = PerformResultAnalysis( resultfilename )
%PerformResultAnalysis 
%Take Tracklet info and analyze for activity measurements

%load result mat file to have variables in function workspace
load(resultfilename);

%only use the completed tracks ( = results for actual fish)
CompletedTracks = tracklets(completeTracks);

%create result struct holding whatever result fields
results.numberoftrackedfish = 0;
results.activitypercentage = [];
results.averagedelay = [];
results.averagedistancetravelled = [];
results.averagedistancetravelledactivefish = [];
try
    %get number of tracked of fish
    results.numberoftrackedfish = numel(completeTracks);
    %calculate activity percentage
    results.activitypercentage = numel(extractfield(CompletedTracks, 'latency'))/size(CompletedTracks,2)*100;
    %calculate average latency
    results.averagedelay = mean(extractfield(CompletedTracks, 'latency'));
    %calculate averagedistancetravelled
    results.averagedistancetravelled = mean(extractfield(CompletedTracks, 'distanceTraveledAfterPulse'));
    %calculate averagedistancetravelled for only active fish
    results.averagedistancetravelledactivefish = mean([CompletedTracks(find(~cellfun(@isempty,{CompletedTracks(:).latency}))).distanceTraveledAfterPulse]);
catch me
    %# report the problematic image, and the reason for failure
    dispstat(me.message, 'keepthis', 'timestamp');
end

%save variables in result file again
save(resultfilename, 'results', '-append');

end

