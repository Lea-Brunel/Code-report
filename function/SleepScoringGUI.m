function SleepScoringGUI
    % Main GUI figure
    hFig = figure('Name', 'Sleep Scoring GUI', 'NumberTitle', 'off', ...
                  'MenuBar', 'none', 'ToolBar', 'none', ...
                  'KeyPressFcn', @keyPressCallback, ...
                  'CloseRequestFcn', @closeGui);

    % Variables
    signalData = [];  % Loaded signal data
    hypnogram = [];   % Hypnogram vector
    time = [];        % Time vector
    currentWindow = [0, 10];  % Current time window in seconds
    samplingRate = 1000;  % Default sampling rate (can be adjusted)
    scoringOnsets = [];  % Stores scoring onsets and their labels
    
    % UI Components
    btnLoad = uicontrol('Style', 'pushbutton', 'String', 'Load LVM File', ...
                        'Units', 'normalized', 'Position', [0.01, 0.95, 0.15, 0.05], ...
                        'Callback', @loadFile);
                    
    btnExport = uicontrol('Style', 'pushbutton', 'String', 'Export Hypnogram', ...
                          'Units', 'normalized', 'Position', [0.17, 0.95, 0.15, 0.05], ...
                          'Callback', @exportHypnogram);

    % Axes for plotting
    ax = axes(hFig, 'Position', [0.1, 0.1, 0.8, 0.8]);
    title(ax, 'Signal Channels'); xlabel(ax, 'Time (s)'); ylabel(ax, 'Amplitude');
    hold(ax, 'on');

    % Load LVM File
    function loadFile(~, ~)
        [fileName, pathName] = uigetfile('*.lvm', 'Select LVM File');
        if isequal(fileName, 0)
            return;
        end
        data = lvm_import(fullfile(pathName, fileName));

        % --------------------------------------------------
        % separate signals into background and signal
        % --------------------------------------------------

        ch = data.Segment1.data(26:end,1);
        sdata = data.Segment1.data(26:end,2);
        
        CA=[];
        q=[];
        demodulated_signal =[];
        
        for k = 1:25
            
            % create logical vectors
            q = (ch == k-1);
            % apply to extract channel
            % REMOVE LINEAR TRENDS
            %srm(:,k) = detrend(sdata(q(:,k)),'linear');
            %s(:,k) = srm(:,k) + 2*min(srm(:,k));
            demodulated_signal(k, :) = sdata(q);
        end
        
        % --- Separate background and Signal in alternating frames ---
        % --- choose odd-even samples
        
        m = 1:numel(demodulated_signal(1,:));
        iseven = rem(m,2) == 0; 
        m_even = m(iseven);
        isodd = rem(m,2) == 1; 
        m_odd = m(isodd);
        
        background  = [demodulated_signal(:,m_odd)];
        signal = demodulated_signal(:,m_even);
        
        if numel(m_even)< numel(m_odd)
            signal = [signal nan(25,1)];
        elseif numel(m_even)> numel(m_odd)
            background = [background nan(25,1)];
        end

        % --------------------------------------------------
        % bin signal to reduce Fs to 5 fps
        % --------------------------------------------------

        % ----------------- 1. NORMALIZE -------------------
        temp = [];
        bin_step = 4;
        samplingRate= 10;
        Fs_bin = samplingRate/bin_step;
        percentile_norm=10;

        temp.s_norm = bin_norm_data(signal', (1:24), bin_step, 1, percentile_norm)';
        temp.b_norm = bin_norm_data(background', (1:24), bin_step, 1, percentile_norm)';
        
        % ----------------- 2. Subtract Baseline -----------
        
        temp.diff = temp.s_norm - temp.b_norm;


        signalData = temp.diff';  % Assuming data structure
        time = (0:(size(signalData, 1) - 1)) / samplingRate;
        hypnogram = zeros(1, length(time));  % Initialize hypnogram vector
        scoringOnsets = [];  % Reset scoring onsets
        currentWindow = [0, 1000];  % Reset window
        plotData();
    end

    % Plot Data
    function plotData
        if isempty(signalData)
            return;
        end
        cla(ax);
        idx = time >= currentWindow(1) & time <= currentWindow(2);
        %plot(ax, time(idx), signalData(idx, 3), 'r', 'DisplayName', 'Channel 1');
        plot(ax, time(idx), signalData(idx, 4), 'g', 'DisplayName', 'Channel 2');
        plot(ax, time(idx), signalData(idx, 5), 'b', 'DisplayName', 'Channel 3');
        legend(ax, 'show');
        xlim(ax, currentWindow);
    end

    % Key Press Callback
    function keyPressCallback(~, event)
        if isempty(signalData)
            return;
        end

        switch event.Key
            case 'leftarrow'
                % Scroll left
                currentWindow = currentWindow - diff(currentWindow) / 2;
                if currentWindow(1) < 0
                    currentWindow = [0, diff(currentWindow)];
                end
                plotData();

            case 'rightarrow'
                % Scroll right
                currentWindow = currentWindow + diff(currentWindow) / 2;
                if currentWindow(2) > max(time)
                    currentWindow = [max(time) - diff(currentWindow), max(time)];
                end
                plotData();

            case {'1', '2', '3'}
                % Sleep scoring
                cp = get(ax, 'CurrentPoint');
                xClick = cp(1, 1);  % x-axis value
                if xClick >= 0 && xClick <= max(time)
                    scoringOnsets = [scoringOnsets; xClick, str2double(event.Key)];
                    idx = round(xClick * samplingRate);
                    hypnogram(idx:end) = str2double(event.Key);
                    disp(['Scored time ', num2str(xClick, '%.2f'), ' as ', event.Key]);
                end
        end
    end

    % Export Hypnogram
    function exportHypnogram(~, ~)
        if isempty(hypnogram)
            msgbox('No hypnogram to export!', 'Error', 'error');
            return;
        end
        [fileName, pathName] = uiputfile('hypnogram.mat', 'Save Hypnogram');
        if isequal(fileName, 0)
            return;
        end
        save(fullfile(pathName, fileName), 'hypnogram', 'scoringOnsets');
        msgbox('Hypnogram exported successfully!', 'Success');
    end

    % Close GUI
    function closeGui(~, ~)
        selection = questdlg('Are you sure you want to close?', ...
                             'Close Request', ...
                             'Yes', 'No', 'No');
        switch selection
            case 'Yes'
                delete(hFig);
            case 'No'
                return;
        end
    end
end
