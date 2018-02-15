folder_name = uigetdir;
searchterm = [folder_name '/*.avi'];
avifiles = dir(searchterm);

for i=1:size(avifiles,1)
    %check for already existing analysis files that do not need to be
    %repeated
    [pathstr,name,ext] = fileparts([folder_name '/' avifiles(i).name]);
    if exist([pathstr '/' name '.mat'], 'file')
        continue;
    end
    %otherwise analyse video
    StartScript([folder_name '/' avifiles(i).name], false);
end

