function formattedNumber = formatWithSpaces(number)
    % Convert the number to a string
    numStr = num2str(number);
    
    % Reverse the string to start inserting spaces from the end
    reverseNumStr = fliplr(numStr);
    
    % Insert spaces every three characters
    parts = regexp(reverseNumStr, '.{1,3}', 'match');
    spacedNumStr = strjoin(parts, ' ');
    
    % Reverse the string back to the original order
    formattedNumber = fliplr(spacedNumStr);
end
