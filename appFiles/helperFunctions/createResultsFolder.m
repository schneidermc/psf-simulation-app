function fullPathSave = createResultsFolder(folderName, parentFolder)
    if nargin < 2
        parentFolder = pwd;
    end
    fullPathSave = fullfile(parentFolder, folderName);
    if ~isfolder(fullPathSave)
        mkdir(fullPathSave)
    end
end