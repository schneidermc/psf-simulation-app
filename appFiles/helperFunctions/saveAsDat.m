function saveAsDat(pathSave,filename,dataMatrix,columnNames)

format = sprintf('%s\n',repmat(' %u',1,size(dataMatrix,2)));
file = fopen([pathSave,'\',filename,'.dat'],'wt');

if nargin == 4
    % Check column number

    if size(columnNames,2)~=size(dataMatrix,2)
        error('Number of column names does not match data matrix!')
    end
    fprintf(file,[strjoin(columnNames,' '),' \n']);
end

fprintf(file,format,dataMatrix.');
fclose(file);

end