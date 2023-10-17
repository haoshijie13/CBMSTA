function sortedVector = replaceNonZeroWithSorted(inputVector)
    % Find the indices and values of non-zero elements
    nonZeroIndices = find(inputVector ~= 0);
    nonZeroValues = inputVector(nonZeroIndices);

    % Sort and remove duplicates from non-zero elements
    sortedUniqueValues = unique(sort(nonZeroValues));

    % Create a mapping table to map original non-zero values to sorted unique values
    valueMap = containers.Map(sortedUniqueValues, 1:length(sortedUniqueValues));

    % Replace non-zero values in the original vector with sorted values based on the mapping table
    for i = 1:length(nonZeroIndices)
        inputVector(nonZeroIndices(i)) = valueMap(nonZeroValues(i));
    end
    
    sortedVector = inputVector;
end

