function loadedMatrix = loadSingleMatrix(filename)
    % Load the .mat file into a structure
    loadedStruct = load(filename);
    
    % Check that the structure only contains one field
    if length(fieldnames(loadedStruct)) == 1
        % Get the fieldname
        field = fieldnames(loadedStruct);

        % Use dynamic fieldnames to assign the field to your chosen variable
        loadedMatrix = loadedStruct.(field{1});
    else
        error('The selected .mat file should contain a single matrix only!')
    end
end
