
paramInputApp()
% Function: paramInputApp
% Purpose: Create a GUI using uifigure to enter multiple analysis parameters
function paramInputApp()

    % Create figure window
    fig = uifigure('Name', 'Enter Analysis Parameters', 'Position', [50 50 1200 800]);

    % Subject ID
    uilabel(fig, 'Text', 'Subject ID:', 'Position', [30 250 100 22]);
    edtSubID = uieditfield(fig, 'text', 'Position', [140 250 200 22], 'Value', 'sub-01');

    % Group (dropdown)
    uilabel(fig, 'Text', 'Group:', 'Position', [30 210 100 22]);
    ddGroup = uidropdown(fig, 'Items', {'control', 'patient'}, 'Position', [140 210 200 22]);

    % Task Name
    uilabel(fig, 'Text', 'Task Name:', 'Position', [30 170 100 22]);
    edtTask = uieditfield(fig, 'text', 'Position', [140 170 200 22], 'Value', 'stroop');

    % Time Window (e.g., [-100 800])
    uilabel(fig, 'Text', 'Time Window:', 'Position', [30 130 100 22]);
    edtTW = uieditfield(fig, 'text', 'Position', [140 130 200 22], 'Value', '[-100 800]');

    % Checkbox: Baseline Correction
    cbBaseline = uicheckbox(fig, 'Text', 'Apply Baseline Correction', 'Position', [140 90 200 22]);

    % Submit button
    btn = uibutton(fig, 'Text', 'Submit', 'Position', [140 40 100 30], 'ButtonPushedFcn', @(btn,event) submitCallback());

    % Callback function
    function submitCallback()
        try
            params.subID = edtSubID.Value;
            params.group = ddGroup.Value;
            params.taskName = edtTask.Value;
            params.timeWindow = str2num(edtTW.Value); %#ok<ST2NM>
            params.doBaseline = cbBaseline.Value;

            % Display to command window or pass to workspace
            assignin('base', 'params', params);
            disp('Parameters saved to workspace variable: params');
            close(fig);
        catch ME
            uialert(fig, ME.message, 'Input Error');
        end
    end
end
