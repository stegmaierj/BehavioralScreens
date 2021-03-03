
%load all result files in folder
%and put them in a new overall result mat file
folder_name = uigetdir;
searchterm = [folder_name '/*.mat'];
if exist([folder_name '/Results.mat'], 'file')
    delete([folder_name '/Results.mat']);
end
matfiles = dir(searchterm);

for i=1:size(matfiles,1)
    [pathstr,name,ext] = fileparts([folder_name '/' matfiles(i).name]);
    load([folder_name '/' matfiles(i).name], 'results');
    Results(i).filename = name;
    tempname = strsplit(name, '_');
    Results(i).time = tempname{2};
    Results(i).name = tempname{3};
    Results(i).numberoftrackedfish = results.numberoftrackedfish;
    Results(i).activitypercentage = results.activitypercentage;
    Results(i).averagedelay = results.averagedelay;
    Results(i).averagedistancetravelled = results.averagedistancetravelled;
    Results(i).averagedistancetravelledactivefish = results.averagedistancetravelledactivefish;
end

%sort for name field according to http://blogs.mathworks.com/pick/2010/09/17/sorting-structure-arrays-based-on-fields/
ResultsFields = fieldnames(Results);
ResultsCell = struct2cell(Results);
sz = size(ResultsCell);
ResultsCell = reshape(ResultsCell, sz(1), []);
ResultsCell = ResultsCell';                       
ResultsCell = sortrows(ResultsCell, 3);
ResultsCell = reshape(ResultsCell', sz);
ResultsSorted = cell2struct(ResultsCell, ResultsFields, 1);

save([pathstr '/' 'Results.mat'],'ResultsSorted');