function [] = fboth(fileHandle, message)

    if strfind(message, 'Processed')
        dispstat(sprintf(message),'timestamp');
    elseif strfind(message, 'Skipped')
        dispstat(sprintf(message),'timestamp');
    elseif strfind(message, 'Processing tracklet')
        dispstat(sprintf(message),'timestamp');
    else
       dispstat(sprintf(message),'keepthis','keepprev','timestamp');
       fprintf(fileHandle, [message '\n']);
    end
    
end