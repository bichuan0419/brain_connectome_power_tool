function iIncrementWaitbar(wb)
    if isempty(wb.UserData)
        wb.UserData = [0, 100]; % Default values, 0 progress out of 100 total
    end
    
    ud = wb.UserData;
    ud(1) = ud(1) + 1; % Increment progress
    
    % Check if the current progress does not exceed total progress
    if ud(1) > ud(2)
        ud(1) = ud(2); % Prevent progress from exceeding total
    end
    
    waitbar(ud(1) / ud(2), wb);
    wb.UserData = ud;
end
