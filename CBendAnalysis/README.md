# CBendAnalysis

# Prerequisites
To be able to run the CBendAnalysis script, make sure to have one of the latest MATLAB versions installed including the Image Processing Toolbox. Note that compressed video files may not be processible if the codices are not installed properly on your system.

# Using the CBendAnalysis Scripts
The main script to start an analysis is named `BatchAnalysis.m`. If you run the script it will ask for an input folder and expects a folder containing `*.avi` video files. All files will then be sequentially processed and the statistics are summarized in a joint `*.csv` file. If you want to process a single video file, you can also simply call the function `AnalyzeSingleVideo(inputFileName, generateResultVideo)` by providing an absolute path to the input video and a boolean flag that indicates wheather to generate an output video or not.

To speedup the processing if only quantitative information is required, you may want to disable the writing of video files in the main script `BatchAnalysis.m` by setting in line 36 to `generateResultVideo = false;` or in case a single video is processed, by setting the second parameter of `AnalyzeSingleVideo(inputFileName, generateResultVideo)` to false.

# Interpreting the Generated Output Tables

For each detected pulse a full list of tracked objects with their properties is generated (`*_PulseXXX_AllTracks.csv`) as well as a filtered one that only includes tracks that span over the entire time interval (`_PulseXXX.csv`). Each of these `*.csv` files contains the following fields:
- Filename: The filename of the input video file.
- Pulse: The index of the pulse (starting at 1 for the first pulse).
- LarvaID: The larva id, which is constant over the entire tracked duration. This ID can also be used to find the corresponding larva in the generated video file.	
- StartTime: The first frame the larva was found and tracked.
- EndTime: The last frame the larva was found. Tracking ends here.	
- NumTrackedFrames: The number of frames the larva was successfully tracked.
- IsComplete: Boolean flag that indicates if the track is only partial (e.g., debris or a partly tracked larva) or if the larva was tracked for the entire duration.
- IsActive: If the larva moves after the pulse was applied it is considered active.
- Latency: The number of frames after the pulse and before the larva starts moving.
- DistanceTraveledBeforePulse: The distance traveled before the pulse measured in pixel.
- DistanceTraveledAfterPulse: The distance traveled after the pulse measured in pixel.
- TotalDistanceTraveled: The total distance traveled in the analyzed video measured in pixel.

In addition to the individual measurements, the `BatchAnalysis.m` script also combines measurements across all analyzed videos. The `combinedResults.csv` sums over all larva and the following features are extracted:
- Filename: The filename of the input video file.
- Pulse: The index of the pulse (starting at 1 for the first pulse).
- NumCompletelyTrackedLarva: The number of larva that were tracked for the entire duration.
- NumActiveLarva: The number of larva that were active during the time interval.
- MeanLatency: The average latency across all analyzed larvae.
- MeanDistanceTraveled: The average distance traveled by all larvae measured in pixel.
- MeanDistanceTraveledActiveLarva: The mean distance traveled by the active larvae measured in pixel.
- SumDistanceTraveledBeforePulse (-FramePadding -> t0): The sum of the distance traveled of all larvae BEFORE the pulse. In the default settings 70 frames before the pulse.
- SumDistanceTraveledAfterPulse (t0 -> FramePadding): The sum of the distance traveled of all larvae AFTER the pulse. In the default settings 70 frames after the pulse.
- ActivityRatio (D_AfterImpulse / D_BeforeImpulse): The ratio of distance traveled after vs. distance traveled before the pulse. A value of 1 indicates no change, values larger than 1, e.g., an n-fold increase of the total movement after the pulse and values smaller than one that movements before the pulse were actually stronger. If there is no movement at all before the pulse, this can result in an Inf/NaN value, due to the division by zero, which is obvious if DistanceTraveledBeforePulse is 0.

 Note that `DistanceTraveledAfterPulse` and `TotalDistanceTraveled` can be different, as the latter includes all frames whereas the former only uses the specified 70 frames after the pulse. This also means that `DistanceTraveledBeforePulse` and `DistanceTraveledAfterPulse` do not necessarily sum to `TotalDistanceTraveled`.
