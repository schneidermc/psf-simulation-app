function [zernikeIndices, zernikeWeights] = loadZernikeIndicesAndWeights(filename)
    [~,~,extension] = fileparts(filename);
    switch extension
        case '.mat'
            % Load the .mat file into a structure
            loadedStruct = load(filename);
            
            % Check that the structure only contains one field
            if length(fieldnames(loadedStruct)) == 1
                % Get the fieldname
                field = fieldnames(loadedStruct);
        
                % Use dynamic fieldnames to assign the field to your chosen variable
                loadedMatrix = loadedStruct.(field{1});
                zernikeIndices = loadedMatrix(:,1);
                zernikeWeights = loadedMatrix(:,2);
            else
                error(['The selected .mat file should contain a single matrix ' ...
                    'containing the Zernike indices as first column and Zernike coefficients as second column!'])
            end
        case '.csv'
            indicesWeights = readmatrix(filename);
            zernikeIndices = indicesWeights(:,1);
            zernikeWeights = indicesWeights(:,2);
        otherwise
            error('Function loadZernikeIndicesAndWeights not implemented for this file extension!')
    end
end
