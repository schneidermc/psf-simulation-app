function s = appendStruct(s1,s2,overwriteDuplicates)
    arguments
        s1 struct
        s2 struct
        overwriteDuplicates logical = 1
    end

    s = [struct2cell(s1); struct2cell(s2)];
    allNames = [fieldnames(s1); fieldnames(s2)];

    if overwriteDuplicates
        [~,indCombined] = unique(allNames,'last');
        s = cell2struct(s(indCombined),allNames(indCombined));
    else
        if length(unique(allNames)) ~= (length(fieldnames(s1)) + length(fieldnames(s2)))
            warning('Duplicated fieldnames!')
        end
        s = cell2struct(s,allNames);
    end
end