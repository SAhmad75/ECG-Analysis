function ECGAnalysisGUIWithDropdownF
    % Specify the path where your ECG data files are located
    dataPath = 'mit-bih-arrhythmia-database-1.0.0/';

    % Get a list of ECG data file names in the specified path
    ecgFiles = dir(fullfile(dataPath, '*.dat'));
    fileNames = {ecgFiles.name};

    % Create main figure
    mainFig = figure('Name', 'ECG Analysis GUI with Dropdown', 'NumberTitle', 'off', 'Position', [100, 100, 600, 500]);
    
    messageLabel = uicontrol('Style', 'text', 'String', 'Select any Signal and perform analysis', ...
        'Position', [150, 420, 350, 20], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');

    % Create dropdown for ECG data files
    fileDropdown = uicontrol('Style', 'popupmenu', 'String',  fileNames, ... ...
        'Position', [50, 370, 200, 30], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.1, 0.9, 0.7]);

    % Create analyze button
    analyzeButton = uicontrol('Style', 'pushbutton', 'String', 'Analyze ECG', ...
        'Position', [50, 330, 120, 30], 'Callback', @analyzeButtonCallback, 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.1, 0.9, 0.7]);

    % Create radio button group
    radioGroup = uibuttongroup('Parent', mainFig, 'Position', [0.05, 0.05, 0.9, 0.6], 'Title', 'Select Image to View', 'FontSize', 12, 'FontWeight', 'bold');

    % Create radio buttons within the group
    radio1 = uicontrol('Style', 'radiobutton', 'Parent', radioGroup, 'String', 'Comparison of QRS Complex Annotations', ...
        'Position', [20, 230, 300, 30], 'FontSize', 10, 'FontWeight', 'bold');
    radio2 = uicontrol('Style', 'radiobutton', 'Parent', radioGroup, 'String', 'ECG with QRS complex and QT interval markings ', ...
        'Position', [20, 200, 300, 30], 'FontSize', 10, 'FontWeight', 'bold');
    radio3 = uicontrol('Style', 'radiobutton', 'Parent', radioGroup, 'String', 'QRS complexes', ...
        'Position', [20, 170, 300, 30], 'FontSize', 10, 'FontWeight', 'bold');
    radio4 = uicontrol('Style', 'radiobutton', 'Parent', radioGroup, 'String', 'ECG with QRS complex markings', ...
        'Position', [20, 140, 300, 30], 'FontSize', 10, 'FontWeight', 'bold');
    radio5 = uicontrol('Style', 'radiobutton', 'Parent', radioGroup, 'String', 'ECG with P-wave markings', ...
        'Position', [20, 110, 300, 30], 'FontSize', 10, 'FontWeight', 'bold');
    radio6 = uicontrol('Style', 'radiobutton', 'Parent', radioGroup, 'String', 'Frequency Domain representation of ECG signal', ...
        'Position', [20, 80, 300, 30], 'FontSize', 10, 'FontWeight', 'bold');
    radio7 = uicontrol('Style', 'radiobutton', 'Parent', radioGroup, 'String', '3D version of signals', ...
        'Position', [20, 50, 300, 30], 'FontSize', 10, 'FontWeight', 'bold');
    radio8 = uicontrol('Style', 'radiobutton', 'Parent', radioGroup, 'String', 'RR intervals for heart rate calculation', ...
        'Position', [20, 20, 300, 30], 'FontSize', 10, 'FontWeight', 'bold');
  
    % Callback function for the analyze button
    function analyzeButtonCallback(~, ~)
        [~, config] = wfdbloadlib;

        % Get the selected file name from the dropdown
        selectedFile = fileNames{get(fileDropdown, 'Value')};
        
        % Load ECG and annotations for the selected file
        N = 10000;
        [ecg, Fs, tm] = rdsamp([dataPath selectedFile], 1, N);

        % Load annotations
        [~, ~, ~, ~, ann] = rdann([dataPath selectedFile(1:end-4)], 'atr', 1, N);
        
        [qrs_amp_raw, qrs_i_raw, ~] = pan_tompkin(ecg, Fs, 0);
        
        [RR, tms]=ann2rr([dataPath selectedFile(1:end-4)], 'atr', 1, N);
        % Call the corresponding analysis function based on the selected radio button
        switch get(get(radioGroup, 'SelectedObject'), 'String')
            case 'Comparison of QRS Complex Annotations'
                analyzeComparison(ecg, tm, ann, qrs_i_raw, qrs_amp_raw, N);
            case 'ECG with QRS complex and QT interval markings '
                analyzeQTInterval(qrs_i_raw, ecg);
            case 'QRS complexes'
                analyzeQRSComplexes(ecg, tm, ann);
            case 'ECG with QRS complex markings'
                analyzeQRSComplexMarkings(ecg, qrs_i_raw, Fs);
            case 'ECG with P-wave markings'
                analyzePWaveMarkings(ecg, qrs_i_raw, Fs);
            case 'Frequency Domain representation of ECG signal'
                analyzeFrequencyDomain(ecg, Fs);
            case '3D version of signals'
                analyze3DVersion(ecg, tms, tm, RR);
            case 'RR intervals for heart rate calculation'
                analyzeRRIntervals(ecg, tm, ann);
        end
    end

    % Analysis functions for each radio button
    function analyzeComparison(ecg, tm, ann, qrs_i_raw, qrs_amp_raw, N)
        % Your existing code for Comparison of QRS Complex Annotations
        figure
        plot(tm(1:N), ecg(1:N)); hold on; grid on
        plot(tm(ann(ann < N) + 1), ecg(ann(ann < N) + 1), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        plot(tm(qrs_i_raw), qrs_amp_raw, 'gx', 'MarkerSize', 8, 'LineWidth', 2);
        legend('ECG Signal', 'Doctor Annotations', 'Pan-Tompkins Algorithm');
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Comparison of QRS Complex Annotations');
        disp('Analysis for Comparison of QRS Complex Annotations');
    end

    function analyzeQTInterval(qrs_i_raw, ecg)
        % Your existing code for ECG with QRS complex and QT interval markings
        % Measure QT interval
        qt_intervals = zeros(1, length(qrs_i_raw));

        % Plot the ECG signal
        figure;
        plot(ecg); hold on;

        for i = 1:length(qrs_i_raw)-1
            % Define a window around the QRS complex
            qrs_start = qrs_i_raw(i);
            qrs_end = qrs_i_raw(i+1);

            % Ensure the window is within the signal boundaries
            qrs_start = max(qrs_start, 1);
            qrs_end = min(qrs_end, length(ecg));

            % Extract the region after the QRS complex (assumption: QT interval follows QRS complex)
            qt_interval = ecg(qrs_i_raw(i)+1 : qrs_end);

            % Identify the end of the T-wave (you may need to adjust this based on your data)
            [~, t_wave_end_index] = min(qt_interval);
            t_wave_end = qrs_i_raw(i) + t_wave_end_index;

            % Calculate QT interval duration
            qt_intervals(i) = t_wave_end - qrs_i_raw(i);

            % Mark QT interval on the ECG signal plot
            plot(qrs_i_raw(i):t_wave_end, ecg(qrs_i_raw(i):t_wave_end), 'g', 'LineWidth', 2);
        end

        % Mark QRS complexes
        scatter(qrs_i_raw, ecg(qrs_i_raw), 'r', 'filled');

        xlabel('Sample Index');
        ylabel('Amplitude');
        title('ECG with QRS Complex and QT Interval Markings');

        % Display QT interval durations
        fprintf('Mean QT Interval: %.2f ms\n', mean(qt_intervals) * 1000);
        fprintf('Standard Deviation of QT Intervals: %.2f ms\n', std(qt_intervals) * 1000);
        disp('Analysis for ECG with QRS complex and QT interval markings');

    end

    function analyzeQRSComplexes(selectedFile)
        % Your existing code for QRS complexes
        % (Omitted for brevity)
        disp('Analysis for QRS complexes');
    end

    function analyzeQRSComplexMarkings(ecg, qrs_i_raw, Fs)
        % Your existing code for ECG with QRS complex markings
        % Assess QRS complex width and morphology
        qrs_widths = zeros(1, length(qrs_i_raw));
        qrs_morphologies = cell(1, length(qrs_i_raw));

        for i = 1:length(qrs_i_raw)
            % Define a window around the QRS complex
            qrs_start = qrs_i_raw(i) - round(0.1*Fs);
            qrs_end = qrs_i_raw(i) + round(0.2*Fs);

            % Ensure the window is within the signal boundaries
            qrs_start = max(qrs_start, 1);
            qrs_end = min(qrs_end, length(ecg));

            % Extract the QRS complex
            qrs_complex = ecg(qrs_start:qrs_end);

            % Calculate QRS complex width
            qrs_widths(i) = length(qrs_complex) / Fs;

            % Store QRS complex morphology
            qrs_morphologies{i} = qrs_complex;
        end

        % Plot the ECG signal with QRS complex markings
        figure;
        plot(ecg); hold on;
        scatter(qrs_i_raw, ecg(qrs_i_raw), 'r', 'filled');
        xlabel('Sample Index');
        ylabel('Amplitude');
        title('ECG with QRS Complex Markings');

        % Display QRS complex widths
        fprintf('Mean QRS Complex Width: %.2f ms\n', mean(qrs_widths) * 1000);
        fprintf('Standard Deviation of QRS Complex Widths: %.2f ms\n', std(qrs_widths) * 1000);

        % Plot a sample of QRS complexes for morphology assessment
        figure;
        for i = 1:min(5, length(qrs_i_raw))
            subplot(5, 1, i);
            plot(qrs_morphologies{i});
            title(['QRS Complex ' num2str(i)]);
        end
        disp('Analysis for ECG with QRS complex markings');
    end

    function analyzePWaveMarkings(ecg, qrs_i_raw)
        % Your existing code for ECG with P-wave markings
        % Extract P-waves and PR intervals
        p_wave_start = zeros(1, length(qrs_i_raw)-1);
        p_wave_end = zeros(1, length(qrs_i_raw)-1);
        pr_intervals = zeros(1, length(qrs_i_raw)-1);

        for i = 1:length(qrs_i_raw)-1
            % Define a window around the QRS complex
            qrs_start = qrs_i_raw(i);
            qrs_end = qrs_i_raw(i+1);

            % Extract the region between QRS complexes
            ecg_segment = ecg(qrs_start:qrs_end);

            % Find the P-wave: you may need to adjust this based on your data
            [~, p_wave_index] = max(ecg_segment(1:round(0.2*Fs)));
            p_wave_start(i) = qrs_start + p_wave_index - 1;

            % Find the end of the P-wave: you may need to adjust this based on your data
            [~, p_wave_end_index] = min(ecg_segment(round(0.2*Fs):end));
            p_wave_end(i) = qrs_start + round(0.2*Fs) + p_wave_end_index - 2;

            % Calculate PR interval
            pr_intervals(i) = p_wave_start(i) - p_wave_end(i);
        end

        % Plot the ECG signal with P-wave markings
        figure;
        plot(ecg); hold on;
        scatter(p_wave_start, ecg(p_wave_start), 'g', 'filled');
        scatter(p_wave_end, ecg(p_wave_end), 'r', 'filled');
        xlabel('Sample Index');
        ylabel('Amplitude');
        title('ECG with P-wave Markings');

        % Display PR intervals
        fprintf('Mean PR Interval: %.2f ms\n', mean(pr_intervals) * 1000);
        fprintf('Standard Deviation of PR Intervals: %.2f ms\n', std(pr_intervals) * 1000);
        disp('Analysis for ECG with P-wave markings');
    end

    function analyzeFrequencyDomain(ecg, Fs)
        % Your existing code for Frequency Domain representation of ECG signal
        % Working in frequency domain

        % Compute the FFT
        NFFT = 2^nextpow2(length(ecg)); % Compute FFT length based on the signal length
        Y = fft(ecg, NFFT); % Compute the FFT of the ECG signal
        Y = Y(1:NFFT/2); % We only need a one-sided FFT plot
        Y_abs = 1/NFFT*abs(Y); % Calculate the magnitude and normalize the spectrum

        % Scale the frequency axis and calculate the corresponding frequencies
        f_fft = (0:NFFT/2-1)*Fs/NFFT;

        % Plot the frequency domain representation
        figure
        plot(f_fft, Y_abs);
        title('Frequency Domain Representation of ECG Signal');
        xlabel('Frequency (Hz)');
        ylabel('Magnitude');

        % Add grid lines for better readability
        grid on;
        disp('Analysis for Frequency Domain representation of ECG signal');
    end

    function analyzeRRIntervals(qrs_i_raw, Fs)
        % Your existing code for RR intervals for heart rate calculation
        % Calculate RR intervals in seconds
        RR_intervals = diff(qrs_i_raw) / Fs;

        % Calculate heart rate in beats per minute (bpm)
        heart_rate = 60 / mean(RR_intervals);

        fprintf('Estimated Heart Rate: %.2f bpm\n', heart_rate);

        % Plot RR intervals for visualization
        figure;
        plot(RR_intervals, '-o');
        xlabel('RR Interval Index');
        ylabel('RR Interval (s)');
        title('RR Intervals for Heart Rate Calculation');
        grid on;
        disp('Analysis for RR intervals for heart rate calculation');
    end

    % Set the callback for radio buttons
    set(radioGroup, 'SelectionChangedFcn', @radioButtonCallback);
end
