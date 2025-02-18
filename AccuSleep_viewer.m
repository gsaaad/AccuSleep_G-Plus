function [message] = AccuSleep_viewer(EEG, EMG, SR, epochLen, userLabels, savepath)
    %AccuSleep_viewer A GUI for manually assigning sleep stage labels to EEG/EMG data.
    %   Zeke Barger 103019
    %   George Saad 021425
    %   Arguments:
    %   EEG: the EEG signal as a vector
    %   EMG: the EMG signal as a vector
    %   SR: the sampling rate for the EEG and EMG in Hz
    %   epochLen: the desired epoch length for sleep stage labels, in seconds.
    %             Values below 2.5 are not supported, and values above 5 are not
    %             recommended.
    %   labels (optional): a vector of sleep stage labels. The length
    %             must match the number of epochs in the EEG/EMG signal.
    %             1 = REM, 2 = wake, 3 = NREM, 4 = undefined
    %   savepath (optional): a filename for saving the sleep stage labels.
    %             AccuSleep_GUI uses this argument, but if you are calling this
    %             function yourself, you can probably ignore it.
    %
    %   Press the 'help' button for the user manual / keyboard shortcuts
    %
    %   Sleep stage labels will be saved to .mat files. The filename can be
    %   anything, but the labels will be stored in a variable called 'labels'
    %
    %   Likewise, labels can be loaded from a .mat file as long as it contains
    %   a variable called 'labels' that is a vector with the same number of
    %   epochs as the EEG/EMG signals and values ranging from 1 to 4.
    
    %% Check the inputs
    G = struct; % holds everything
    
    % make sure we have at least 3 arguments
    switch nargin
        case {0,1,2,3}
            disp('Not enough arguments')
            message = 'ERROR: Not enough arguments';
            return
        case 4
            G.labels = [];
        case {5,6}
            if ~isempty(userLabels) % user tried to pass some labels
                % check whether they have the correct format
                q = struct;
                q.labels = userLabels;
                if checkLabels(q, 0) % meets certain conditions
                    G.labels = q.labels;
                    if isrow(G.labels)
                        G.labels = G.labels';
                    end
                else
                    message = 'ERROR: Labels file is incorrectly formatted';
                    return
                end
            else
                G.labels = [];
            end
    end
    
    G.originalSR = SR; % EEG/EMG sampling rate
    G.SR = 256; % sampling rate used when calculating spectrogram and processed EMG
    G.epochLen  = epochLen; % length of one epoch (spectrogram column) in seconds
    
    if length(EEG) ~= length(EMG)
        message = 'ERROR: EEG and EMG are different lengths';
        return
    end
    
    % Filter the EEG and EMG data 
    EEG = lowpass(EEG, 25, SR); 
    %EEG = bandpass(EEG,[.3 25], SR);
    EMG = bandpass(EMG,[25 50], SR);
    
    G.EEG = EEG - mean(EEG);
    G.EMG = EMG - mean(EMG);
    clear('EMG','EEG');
    
    % create spectrogram and process EMG at a standard SR (128)
    [spec, tAxis, fAxis] = createSpectrogram(standardizeSR(G.EEG, G.originalSR, G.SR), G.SR, G.epochLen);
    G.processedEMG = processEMG(standardizeSR(G.EMG, G.originalSR, G.SR), G.SR, G.epochLen);
    % set ceiling for EMG trace at 2.5 SD when plotting
    G.cappedEMG = G.processedEMG;
    emgCap = mean(G.cappedEMG) + 2.5*std(G.cappedEMG);
    G.cappedEMG(G.cappedEMG > emgCap) = emgCap;
    
    % set various parameters
    G.show = 15; %  default number of bins to display on screen
    G.dt = 1/G.originalSR; % duration of each EEG/EMG sample in seconds
    G.advance = 0; % whether to advance automatically when a state is assigned
    G.colors = [1 1 1; .47 .67 .19; .14 .2 .57; 0.996 0.758 0.039; 0 0 0]; %colors for sleep stages
    G.mid = ceil(G.show/2); % important for plotting the current time marker - middle of G.show
    G.savepath = ''; % where to save the sleep stage labels
    G.nbins = length(tAxis); % total number of time bins in the recording
    G.unsavedChanges = 0; % whether anything has changed since user last saved
    
    % check to make sure labels has the proper length
    if ~isempty(G.labels)
        if length(G.labels) ~= G.nbins
            disp(['length of labels must be ',num2str(G.nbins),' for this recording'])
            message = 'ERROR: Length of labels file does not match EEG/EMG. Check SR / epoch length?';
            return
        end
    else % if no sleep stages provided, set to undefined
        G.labels = ones(G.nbins,1) * 4;
    end
    
    % get spectrogram and time axes
    showFreqs = find(fAxis <= 30); % only show frequencies under 30 Hz
    G.specTs = (1:G.nbins)*G.epochLen - G.epochLen/2; % spectrogram time axis, in seconds
    G.specTh = G.specTs./3600; % spectrogram time axis, in hours
    G.spectrogram = spec(:,showFreqs); % our EEG spectrogram
    if nargin == 6 % this probably only happens when called by AccuSleep_GUI
        G.f = fAxis; % frequency axis
        G.s = spec;
        G.savepath = savepath;
    end
    
    % take a sample of the spectrogram to help initialize the colormap
    sampleBins = randperm(G.nbins, round(G.nbins/10));
    specSample = reshape(spec(sampleBins,showFreqs),1,length(sampleBins)*length(showFreqs));
    G.caxis1 = prctile(specSample,[6 98]);
    G.cmax = G.caxis1(2);
    clear spec
    
    % figure out how to scale eeg and emg visually
    G.eegLen = length(G.EEG); % length of our recording, in samples
    yl = prctile(abs(G.EEG(1:10:end)),95);
    G.eegYlim = [-2.2*yl, 2.2*yl];
    G.emgYlim = G.eegYlim;
    
    % load our colormap
    try
        G.colormap = AccuSleep_colormap();
    catch
        try
            G.colormap = parula;
        catch
            G.colormap = jet;
        end
    end
    
    
    
    %% Make the figure window
    
    WIN = figure('Units', 'Normalized', 'CloseRequestFcn',@closeReq,...
        'Position', [0.08, 0.12, 0.83, 0.75],'KeyPressFcn',@keypress,...
        'Menubar', 'none','Color', 'w', 'Name', 'AccuSleep_viewer');
    
    % create axes
    % panel divider and y labels
    G.A5 = axes('Units', 'Normalized', 'Position', [0 .579 .05 .02],'XColor','w','YColor','w');
    hold(G.A5,'on');
    set(G.A5,'Box','off','XLim',[0 1],'YLim',[0 0.1],'XTick',[],'YTick',[],'Clipping','off')
    G.A5.Toolbar.Visible = 'off';
    
    % EEG spectrogram
    G.A3 = axes('Units', 'Normalized', 'Position', [0.05 0.75 0.87 0.11]);
    G.A3.Toolbar.Visible = 'off';
    ylabel(G.A3, 'Spec.');
    set(gca, 'FontSize', 10, 'LineWidth', 2, 'XTick', [], 'YTick', []);
    
    % Log Emg
    G.A4 = axes('Units', 'Normalized', 'Position', [0.05 0.66  0.87 0.07]);
    G.A4.Toolbar.Visible = 'off';
    ylabel(G.A4, 'Filtered EMG(20-50)');
    axis(G.A4, 'off');
    
    % EEG signal - 1 Minute
    G.A6a = axes('Units', 'Normalized', 'Position', [0.05 0.52, 0.87 .11]);
    G.A6a.Toolbar.Visible = 'off';
    set(G.A6a,'XTick',[])
    ylabel(G.A6a, 'EEG1');
    
    % % EMG signal - 1 Minute
    G.A7a = axes('Units', 'Normalized', 'Position', [0.05 0.35 0.87 .11]);
    G.A7a.Toolbar.Visible = 'off';
    ylabel(G.A7a, 'EMG1');
    
    % EEG signal
    G.A6 = axes('Units', 'Normalized', 'Position', [0.05 0.18 0.87 .11]);
    G.A6.Toolbar.Visible = 'off';
    ylabel(G.A6, 'EEG2');
    
    % EMG signal
    G.A7 = axes('Units', 'Normalized', 'Position', [0.05 0.03 0.87 .11]);
    G.A7.Toolbar.Visible = 'off';
    ylabel(G.A7, 'EMG2');
    
    % Upper sleep stage labels
    G.A1 = axes('Units', 'Normalized', 'Position', [0.05 0.915 0.87 0.08]);
    G.A1.Toolbar.Visible = 'off';
    % Rotated label
    ylabel(G.A1, 'State');
    
    % Time point indicator
    G.A2 = axes('Units', 'Normalized', 'Position', [0.05 0.895  0.87 0.02],'XTick',[],'YTick',[]);
    G.A2.Toolbar.Visible = 'off';
    
    
    % axis labels - old
    % text(G.A5,0.05,1.87,'State')
    % text(G.A5,0.05,-.1,'EEG1(mV)')
    % text(G.A5,0.05,-.9,'EMG1(mV)')
    % text(G.A5,0.05,-1.6,'EEG2(mV)')
    % text(G.A5,0.05,-2.35,'EMG2(mV)')
    
    linkaxes([G.A1, G.A2, G.A3, G.A4], 'x'); % upper panel x axes should stay linked
    
    % buttons
    G.helpbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized','BackgroundColor',[1 .8 .8],...
        'Position',[.93 .94 .062 .055],'String','Help','Callback',@showHelp,'FontSize',9,...
        'ToolTip','Show help menu (H)');
    G.savebtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized','BackgroundColor',[.8 1 .8],...
        'Position',[.93 .885 .062 .045],'String','Save labels','Callback',@saveCallback,'FontSize',9,...
        'ToolTip','Save labels to file (F)');
    G.loadbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
        'Position',[.93 .845 .062 .03],'String','Load labels','Callback',@loadFile,'FontSize',9,...
        'ToolTip','Load labels from file');
    G.brightbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
        'Position',[.93 .8 .062 .025],'String','Brighter','Callback',@brightSpect,...
        'FontSize',9,'ToolTip','Make EEG spectrogram brighter');
    G.dimbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
        'Position',[.93 .77 .062 .025],'String','Dimmer','Callback',@dimSpect,...
        'FontSize',9,'ToolTip','Make EEG spectrogram dimmer');
    G.selectbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
        'Position',[.93 .705 .062 .04],'String','<html>Select<br>timepoint',...
        'Callback',@fct_selectlocation,'FontSize',9,'ToolTip','Click a timepoint to show (A)');
    G.zoominbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
        'Position',[.93 .66 .062 .025],'String','Zoom IN','Callback',@fct_zoomin_t,...
        'FontSize',9,'ToolTip','Increase zoom level (+)');
    G.zoomoutbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
        'Position',[.93 .63 .062 .025],'String','Zoom OUT','Callback',@fct_zoomout_t,...
        'FontSize',9,'ToolTip','Decrease zoom level (-)');
    % % Create a text box for EEG amplitude display, placed above or below the buttons:
    % G.EEGprev_ampDisplay = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
    %     'Position', [.925 .4 .07 .035], 'String', '-- to --', ...
    %     'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
    % % Create a text box for EEG amplitude display, placed above or below the buttons:
    % G.EEG_ampDisplay = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
    %     'Position', [.925 .295 .07 .035], 'String', '-- to --', ...
    %     'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
    % % Create a text box for EMG amplitude display, placed above or below the buttons:
    % G.EMG_ampDisplay = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
    %     'Position', [.925 .18 .07 .035], 'String', '-- to --', ...
    %     'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
    G.zoomresetbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
        'Position',[.93 .595 .062 .025],'String','Reset zoom','Callback',@fct_zoomreset_t,...
        'FontSize',9,'ToolTip','Reset zoom level (0)');
    G.gui_zoominEEG = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
        'Position',[0.93 0.505 0.02 0.02],'Callback', @fct_zoominEEG,'String','+',...
        'FontSize',9);
    G.gui_zoomoutEEG = uicontrol(WIN,'Style','pushbutton','Units','normalized',...
        'Position',[0.93 0.475 0.02 0.02],'Callback', @fct_zoomoutEEG,'String','-',...
        'FontSize',9);
    G.gui_shiftupEEG = uicontrol(WIN,'Style','pushbutton','Units','normalized',...
        'Position',[0.96 0.475 0.02 0.02],'Callback', @fct_shiftupEEG,'String','\/',...
        'FontSize',9);
    G.gui_shiftdownEEG = uicontrol(WIN,'Style','pushbutton','Units','normalized',...
        'Position',[0.96 0.505 0.02 0.02],'Callback', @fct_shiftdownEEG,'String','/\',...
        'FontSize',9);
    G.gui_zoominEMG = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
        'Position',[0.93 0.28 0.02 0.02],'Callback', @fct_zoominEMG,'String','+',...
        'FontSize',9);
    G.gui_zoomoutEMG = uicontrol(WIN,'Style','pushbutton','Units','normalized',...
        'Position',[0.93 0.25 0.02 0.02],'Callback', @fct_zoomoutEMG,'String','-',...
        'FontSize',9);
    G.gui_shiftupEMG = uicontrol(WIN,'Style','pushbutton','Units','normalized',...
        'Position',[0.96 0.25 0.02 0.02],'Callback', @fct_shiftupEMG,'String','\/',...
        'FontSize',9);
    G.gui_shiftdownEMG = uicontrol(WIN,'Style','pushbutton','Units','normalized',...
        'Position',[0.96 0.28 0.02 0.02],'Callback', @fct_shiftdownEMG,'String','/\',...
        'FontSize',9);
    G.showMenu = uicontrol(WIN,'Style','popupmenu','Units','normalized',...
        'Position',[0.93 0.54 0.062 0.02],'Callback', @fct_showmenu,...
        'String',{'Show 1 epoch','Show 3 epochs','Show 5 epochs','Show 7 epochs','Show 9 epochs','Show 15 epochs'},...
        'Value',3);
    G.rangebtn = uicontrol(WIN,'Style','pushbutton','Units','normalized','BackgroundColor',[.92 .92 .92],...
        'Position',[.93 .16 .062 .025],'String','set range','Callback',@setRange,...
        'FontSize',9,'ToolTip',sprintf(['Set state for range of timepoints (*)',...
        '\nDraw an ROI on the upper sleep stage panel,\nand double-click it']));
    G.nrembtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
        'Position',[.93 .13 .062 .025],'String','NREM','Callback',@(src,evnt)setState(src,evnt,3),...
        'FontSize',9,'ToolTip','Set state to NREM (S)','BackgroundColor',[1 .96 .82]); %.8 .8 .8
    G.wakebtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
        'Position',[.93 .1 .062 .025],'String','Wake','Callback',@(src,evnt)setState(src,evnt,2),...
        'FontSize',9,'ToolTip','Set state to wake (W)','BackgroundColor',[.86 .88 1]); % .93 .5 .93
    G.rembtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
        'Position',[.93 .07 .062 .025],'String','REM','Callback',@(src,evnt)setState(src,evnt,1),...
        'FontSize',9,'ToolTip','Set state to REM (R)','BackgroundColor',[.84 .92 .73]); % .5 1 1
    G.autobox = uicontrol(WIN,'Style','checkbox', 'Units','normalized',...
        'Position',[.93 .005 .062 .06],'String','<html>Auto-<br>scroll','Callback',@scrollCallback,...
        'FontSize',9,'ToolTip','Advance to next time step after assigning label (insert)');
    
    
    % keep track of the current timepoint
    G.index = 1; % index of current time point
    G.timepointS = G.specTs(G.index); % current time point, in seconds
    G.timepointH = G.specTh(G.index); % current time point, in hours
    
    % plot spectrogram
    caxis(G.A3,G.caxis1);
    imagesc(G.A3,G.specTh, fAxis(showFreqs), G.spectrogram',G.caxis1);
    axis(G.A3, 'xy')
    colormap(G.A3,G.colormap);
    G.lims = xlim(G.A3); % store maximum x limits for the upper panel plots
    set(G.A3, 'YTick', 0:5:30, 'YTickLabel', 0:5:30);
    
    % plot processed EMG
    plot(G.A4,G.specTh,G.cappedEMG,'k')
    yr = max(G.cappedEMG) - min(G.cappedEMG); % adjust y limits
    set(G.A4,'XTick',[],'YTick',[],'box','off',...
        'YLim',[min(G.cappedEMG) - .02*yr, max(G.cappedEMG) + .02*yr])
    % y_min = min(G.cappedEMG) - 0.02 * yr;        % add 2% margin below
    % y_max = max(G.cappedEMG) + 0.02 * yr;        % add 2% margin above
    
    % % Optionally, round the tick values to one decimal place:
    % customTicks = round(linspace(y_min, y_max, 6), 1);
    % 
    % % Set the y-axis properties on the axes handle G.A4:
    % set(G.A4, 'YLim', [y_min, y_max], 'YTick', customTicks, 'YTickMode', 'manual');
    % 
    % % Optionally, display these values for verification:
    % fprintf('Dynamic y-axis limits: [%.2f, %.2f]\n', y_min, y_max);
    % fprintf('Custom YTick values: %s\n', mat2str(customTicks));
    
    % Upper sleep stages
    box(G.A1, 'on');
    xlim(G.A1,[G.specTh(1)-G.epochLen/3600, G.specTh(end)-G.epochLen/3600]);
    updateState;
    
    % Plot everything else
    updatePlots;
    axes(G.A1);
    
    message = 'Data loaded successfully';
    
    %% Functions used by buttons, keypresses, etc.
        function updatePlots(~, ~) % update plots when something changes
            
            
            % plot sleep stage in the lower panel
            n = (G.show-1)/2; % number of bins on either side of center to show
            tp = G.timepointS; % time in seconds at the center of the screen
            gi = G.index; % index of bin in the center of the screen
            if G.index < G.mid
                gi = G.mid;
                tp = gi*G.epochLen-G.epochLen/2;
            end
            if G.nbins - G.index < (G.mid-1)
                gi = G.nbins-(G.mid-1);
                tp = gi*G.epochLen-G.epochLen/2;
            end
    
            yr = max(G.cappedEMG) - min(G.cappedEMG);  % range of the data
        
            
            % plot EEG and EMG
            n = round((G.show*G.epochLen)/G.dt/2); % number of samples on either side to show
            i = round(tp / G.dt);
            ii = i-n:i+n; % choose indices to show
            t = tp-n*G.dt:G.dt:tp+n*G.dt;
            ii(ii<=0) = 1;
            ii(ii>=G.eegLen) = G.eegLen;
            
            % --- Plot Processed EMG on G.A7 with Dynamic Y-Axis ---
            cla(G.A7);
            hold(G.A7, 'on');
            xlim(G.A7, [t(1)-G.dt, t(end)]);
            
            % Compute dynamic y-axis limits from the current EMG window data (in volts)
            curEMG = G.EMG(ii);
            y_min_EMG = min(curEMG);
            y_max_EMG = max(curEMG);
            yr_EMG = y_max_EMG - y_min_EMG;
            padding_EMG = 0.02 * yr_EMG;
            
            if yr_EMG < eps
                y_min_dyn_EMG = y_min_EMG - 0.1;
                y_max_dyn_EMG = y_max_EMG + 0.1;
            else
                y_min_dyn_EMG = y_min_EMG - padding_EMG;
                y_max_dyn_EMG = y_max_EMG + padding_EMG;
            end
            
            if y_min_dyn_EMG >= y_max_dyn_EMG
                y_max_dyn_EMG = y_min_dyn_EMG + 0.1;
            end
            
            % Apply the dynamic y-axis limits (data remain in volts)
            % ylim(G.A7, [y_min_dyn_EMG, y_max_dyn_EMG]);
            % set(G.A7, 'XLimMode','manual', 'YLimMode','manual');
    
    
            
            % Plot the EMG data (in volts)
            globalMin_EMG = min(G.EMG);  
            globalMax_EMG = max(G.EMG);  
            globalRange_EMG = globalMax_EMG - globalMin_EMG;
            padding_EMG = 0.02 * globalRange_EMG;  
            
            stableMin_EMG = globalMin_EMG - padding_EMG;
            stableMax_EMG = globalMax_EMG + padding_EMG;
            
            % --- In updatePlots, apply the stable EEG axis limits ---
            set(G.A7, 'YLim', [stableMin_EMG, stableMax_EMG], 'YLimMode', 'manual');
            line(G.A7, t, curEMG, 'Color', 'k', 'LineWidth', 1);
            % Plot red indicator lines for current epoch boundaries:
            indicatorHeight_EMG = 0.1 * (globalRange_EMG);  % Height based on current axis range
            
            % Vertical line at left epoch boundary (starts at bottom, points upward)
            line(G.A7, ones(1,2)*(G.timepointS - G.epochLen/2),[globalMin_EMG, globalMin_EMG + indicatorHeight_EMG],'Color','r', 'LineWidth', 0.5);
            
            % Vertical line at right epoch boundary (starts at bottom, points upward)
            line(G.A7, ones(1,2)*(G.timepointS + G.epochLen/2),[globalMin_EMG, globalMin_EMG + indicatorHeight_EMG],'Color','r', 'LineWidth', 0.5);
            
    
            % Horizontal line connecting the tops of the vertical markers
            line(G.A7, [G.timepointS - G.epochLen/2, G.timepointS + G.epochLen/2], [globalMin_EMG, globalMin_EMG], 'Color','r', 'LineWidth', 0.5);
    
    
        
    
            % 
            % desiredNumTicks = 4;
            % precision = 2;  % display precision: one decimal (in mV)
            % 
            % % Convert dynamic limits from volts to mV:
            % y_min_dyn_EMG_mV = y_min_dyn_EMG * 1000;
            % y_max_dyn_EMG_mV = y_max_dyn_EMG * 1000;
            % 
            % % Generate 4 evenly spaced tick marks between dynamic limits (in mV)
            % customTicks_mV = linspace(y_min_dyn_EMG_mV, y_max_dyn_EMG_mV, desiredNumTicks);
            % 
            % % Force endpoints to match the dynamic limits exactly (in mV)
            % customTicks_mV(1) = round(y_min_dyn_EMG_mV, precision);
            % customTicks_mV(end) = round(y_max_dyn_EMG_mV, precision);
            % 
            % % Ensure final tick values are rounded to one decimal (in mV)
            % customTicks_mV = round(customTicks_mV, precision);
            % EMG_amp_min_mV = customTicks_mV(1) * 1000;
            % EMG_amp_max_mV = customTicks_mV(end) * 1000;
            % EMG_amplitudeStr = sprintf('%d to %d mV', EMG_amp_min_mV, EMG_amp_max_mV);
            % % set(G.EMG_ampDisplay, 'String', EMG_amplitudeStr);
            % 
            % % Check if the tick vector is strictly increasing and has no duplicates
            % if numel(unique(customTicks_mV)) < desiredNumTicks || any(diff(customTicks_mV) <= 0)
            %     % Fallback to predefined ticks that are valid after rounding
            %     customTicks_mV = [-0.2, -0.1, 0.1, 0.2];
            %     customTicks_mV = round(customTicks_mV, precision);
            % end
            % 
            % % Set the YTick values on G.A7 (convert mV back to volts)
            % set(G.A7, 'YTick', customTicks_mV / 1000, 'YTickMode', 'manual');
            % 
            % % Format tick labels for display (in mV)
            % formattedTicks_EMG = arrayfun(@(x) sprintf('%.1f', x), customTicks_mV, 'UniformOutput', false);
            % set(G.A7, 'YTickLabel', formattedTicks_EMG, 'YTickMode', 'manual');
            % 
            % % Verification output:
            % currentYLim_EMG = get(G.A7, 'YLim');
            % disp(['Dynamic EMG Y-axis limits: ', num2str(currentYLim_EMG)]);
            % disp(['Underlying EMG YTick values: ', mat2str(customTicks_mV)]);
            % disp(['Displayed EMG YTick labels: ', strjoin(formattedTicks_EMG, ', ')]);
    
    
    
            % --- Compute Global EMG Limits (in volts) ---
            globalMin_EEG = min(G.EEG);  
            globalMax_EEG = max(G.EEG);  
            globalRange_EEG = globalMax_EEG - globalMin_EEG;
            padding_EEG = 0.02 * globalRange_EEG;  
            
            stableMin_EEG = globalMin_EEG - padding_EEG;
            stableMax_EEG = globalMax_EEG + padding_EEG;
            
            % --- In updatePlots, apply the stable EEG axis limits ---
            set(G.A6, 'YLim', [stableMin_EEG, stableMax_EEG], 'YLimMode', 'manual');
            set(G.A6,'XTick',[])
            
            
            cla(G.A6)
            hold(G.A6, 'on');
            xlim(G.A6,[t(1)-G.dt t(end)]);
            % Compute dynamic y-axis limits based on current EEG window data (in volts)
    
            curEEG = G.EEG(ii);
            y_min_EEG = min(curEEG);
            y_max_EEG = max(curEEG);
            yr_EEG = y_max_EEG - y_min_EEG;
            padding_EEG = 0.02 * yr_EEG;
            if yr_EEG < eps
                y_min_dyn_EEG = y_min_EEG - 0.1;
                y_max_dyn_EEG = y_max_EEG + 0.1;
            else
                y_min_dyn_EEG = y_min_EEG - padding_EEG;
                y_max_dyn_EEG = y_max_EEG + padding_EEG;
            end
            if y_min_dyn_EEG >= y_max_dyn_EEG
                y_max_dyn_EEG = y_min_dyn_EEG + 0.1;
            end
            % Apply the dynamic y-axis limits (still in volts)
            % Plot the EEG data (in volts)
            line(G.A6, t, curEEG, 'Color','k', 'LineWidth', 1);
            indicatorHeight_EEG = 0.1 * (globalRange_EEG);
            
            % Draw vertical red lines at the left and right boundaries of the current epoch,
            % with the line extending from the top (y_max_dyn_EEG) down by indicatorHeight_EEG.
            line(G.A6, ones(1,2) * (G.timepointS - G.epochLen/2), [globalMax_EEG, globalMax_EEG - indicatorHeight_EEG], 'Color','r', 'LineWidth', 0.5);
            line(G.A6, ones(1,2) * (G.timepointS + G.epochLen/2), [globalMax_EEG, globalMax_EEG - indicatorHeight_EEG], 'Color','r', 'LineWidth', 0.5);
            % Draw a horizontal red line across the top of the current epoch:
            line(G.A6, [G.timepointS - G.epochLen/2, G.timepointS + G.epochLen/2], [globalMax_EEG, globalMax_EEG], 'Color','r', 'LineWidth', 0.5);
    
    
            desiredNumTicks = 4;
            tickPrecision = 1;
            y_min_dyn_EEG_mV = y_min_dyn_EEG * 1000;
            y_max_dyn_EEG_mV = y_max_dyn_EEG * 1000;
            customTicks_EEG_mV = linspace(y_min_dyn_EEG_mV, y_max_dyn_EEG_mV, desiredNumTicks);
            customTicks_EEG_mV(1) = round(y_min_dyn_EEG_mV, tickPrecision);
            customTicks_EEG_mV(end) = round(y_max_dyn_EEG_mV, tickPrecision);
            EEG_amp_min_mV = customTicks_EEG_mV(1) * 1000;
            EEG_amp_max_mV = customTicks_EEG_mV(end) * 1000;
            EEG_amplitudeStr = sprintf('%d to %d mV', EEG_amp_min_mV, EEG_amp_max_mV);
            % set(G.EEG_ampDisplay, 'String', EEG_amplitudeStr);
    
            if numel(unique(customTicks_EEG_mV)) < desiredNumTicks || any(diff(customTicks_EEG_mV) <= 0)
                delta = (y_max_dyn_EEG_mV - y_min_dyn_EEG_mV) / (desiredNumTicks - 1);
                customTicks_EEG_mV = y_min_dyn_EEG_mV : delta : y_max_dyn_EEG_mV;
                if numel(customTicks_EEG_mV) ~= desiredNumTicks
                    customTicks_EEG_mV = linspace(y_min_dyn_EEG_mV, y_max_dyn_EEG_mV, desiredNumTicks);
                end
            end
    
    
            % Set the YTick values on G.A6.
            % Since the data are in volts, convert the mV tick values back to volts:
            % set(G.A6, 'YTick', customTicks_EEG_mV / 1000, 'YTickMode', 'manual');
            % % For display, create formatted tick labels (in mV)
            % formattedTicks_EEG = arrayfun(@(x) sprintf('%.1f', x), customTicks_EEG_mV, 'UniformOutput', false);
            % set(G.A6, 'YTickLabel', formattedTicks_EEG, 'YTickMode', 'manual');
            % % For verification, display the dynamic limits and tick values:
            % currentYLim_EEG = get(G.A6, 'YLim');
            % disp(['Dynamic EEG Y-axis limits: ', num2str(currentYLim_EEG)]);
            % disp(['Underlying EEG YTick values: ', mat2str(customTicks_EEG_mV)]);
            % disp(['Displayed EEG YTick labels: ', strjoin(formattedTicks_EEG, ', ')]);
            % Apply the dynamic y-axis limits (data remain in volts)
    
              
    
    
    
                    %% Plot EEG signal - Previous Minute in g.A6a
            % Determine if there is a previous minute available:
            %% Plot EEG signal - Previous Minute in g.A6a (with dynamic YTick generation)
            % --- Plot EEG Signal - Previous Minute in g.A6a (with dynamic YTick generation) ---
            if G.timepointS > 60
                tp_prev = G.timepointS - 60;  % center time for previous minute
                n_prev = round((G.show * G.epochLen) / G.dt / 2);
                i_center_prev = round(tp_prev / G.dt);
                ii_prev = i_center_prev - n_prev : i_center_prev + n_prev;
                ii_prev(ii_prev <= 0) = 1;
                ii_prev(ii_prev >= G.eegLen) = G.eegLen;
                t_prev = (tp_prev - n_prev * G.dt) : G.dt : (tp_prev + n_prev * G.dt);
                
                % EEG-previous minutes -1 minute
                cla(G.A6a);
                hold(G.A6a, 'on');
                xlim(G.A6a, [t_prev(1)-G.dt, t_prev(end)]);
                
                % Compute dynamic y-axis limits for previous EEG (in volts)
                curEEG_prev = G.EEG(ii_prev);
                y_min_EEG_prev = min(curEEG_prev);
                y_max_EEG_prev = max(curEEG_prev);
                yr_EEG_prev = y_max_EEG_prev - y_min_EEG_prev;
                padding_EEG_prev = 0.02 * yr_EEG_prev;
                if yr_EEG_prev < eps
                    y_min_dyn_EEG_prev = y_min_EEG_prev - 0.1;
                    y_max_dyn_EEG_prev = y_max_EEG_prev + 0.1;
                else
                    y_min_dyn_EEG_prev = y_min_EEG_prev - padding_EEG_prev;
                    y_max_dyn_EEG_prev = y_max_EEG_prev + padding_EEG_prev;
                end
                if y_min_dyn_EEG_prev >= y_max_dyn_EEG_prev
                    y_max_dyn_EEG_prev = y_min_dyn_EEG_prev + 0.1;
                end
                ylim(G.A6a, [y_min_dyn_EEG_prev, y_max_dyn_EEG_prev]);
                set(G.A6a, 'XLimMode','manual', 'YLimMode','manual');
                
                % Plot the previous EEG data (in volts)
                line(G.A6a, t_prev, curEEG_prev, 'Color','k', 'LineWidth', 1);
                
                % --- In updatePlots, apply the stable EEG axis limits ---
                set(G.A6a, 'YLim', [stableMin_EEG, stableMax_EEG], 'YLimMode', 'manual');
                
                % --- Dynamically Set YTick Values for g.A6a ---
                % desiredNumTicks = 4;
                % tickPrecision = 1;
                % y_min_dyn_EEG_prev_mV = y_min_dyn_EEG_prev * 1000;
                % y_max_dyn_EEG_prev_mV = y_max_dyn_EEG_prev * 1000;
                % customTicks_EEG_prev_mV = linspace(y_min_dyn_EEG_prev_mV, y_max_dyn_EEG_prev_mV, desiredNumTicks);
                % customTicks_EEG_prev_mV(1) = round(y_min_dyn_EEG_prev_mV, tickPrecision);
                % customTicks_EEG_prev_mV(end) = round(y_max_dyn_EEG_prev_mV, tickPrecision);
                % EEGprev_amp_min_mV = customTicks_EEG_prev_mV(1) * 1000;
                % EEGprev_amp_max_mV = customTicks_EEG_prev_mV(end) * 1000;
                % EEGprev_amplitudeStr = sprintf('%d to %d mV', EEGprev_amp_min_mV, EEGprev_amp_max_mV);
                % % set(G.EEGprev_ampDisplay, 'String', EEGprev_amplitudeStr);
                % if numel(unique(customTicks_EEG_prev_mV)) < desiredNumTicks || any(diff(customTicks_EEG_prev_mV) <= 0)
                %     delta = (y_max_dyn_EEG_prev_mV - y_min_dyn_EEG_prev_mV) / (desiredNumTicks - 1);
                %     customTicks_EEG_prev_mV = y_min_dyn_EEG_prev_mV : delta : y_max_dyn_EEG_prev_mV;
                %     if numel(customTicks_EEG_prev_mV) ~= desiredNumTicks
                %         customTicks_EEG_prev_mV = linspace(y_min_dyn_EEG_prev_mV, y_max_dyn_EEG_prev_mV, desiredNumTicks);
                %     end
                % end
                % customTicks_EEG_prev_mV = round(customTicks_EEG_prev_mV, tickPrecision);
                % set(g.A6a, 'YTick', customTicks_EEG_prev_mV / 1000, 'YTickMode', 'manual');
                % % For display, create formatted tick labels (in mV)
                % formattedTicks_EEG_prev = arrayfun(@(x) sprintf('%.1f', x), customTicks_EEG_prev_mV, 'UniformOutput', false);
                % set(g.A6a, 'YTickLabel', formattedTicks_EEG_prev, 'YTickMode', 'manual');
                % set(g.A6a,'XTick',[])
    
                % Verification output for previous minute:
                % currentYLim_EEG_prev = get(g.A6a, 'YLim');
                % disp(['Dynamic EEG Y-axis limits for previous minute: ', num2str(currentYLim_EEG_prev)]);
                % disp(['Underlying EEG YTick values for previous minute: ', mat2str(customTicks_EEG_prev_mV)]);
                % disp(['Displayed EEG YTick labels for previous minute: ', strjoin(formattedTicks_EEG_prev, ', ')]);
    
    
                % EMG-previous minutes -1 minute
                cla(G.A7a);
                hold(G.A7a, 'on');
                xlim(G.A7a, [t_prev(1)-G.dt, t_prev(end)]);
                
                % Compute dynamic y-axis limits for previous EEG (in volts)
                curEMG_prev = G.EMG(ii_prev);
                y_min_EMG_prev = min(curEMG_prev);
                y_max_EMG_prev = max(curEMG_prev);
                yr_EMG_prev = y_max_EMG_prev - y_min_EMG_prev;
                padding_EMG_prev = 0.02 * yr_EMG_prev;
                if yr_EMG_prev < eps
                    y_min_dyn_EMG_prev = y_min_EMG_prev - 0.1;
                    y_max_dyn_EMG_prev = y_max_EMG_prev + 0.1;
                else
                    y_min_dyn_EMG_prev = y_min_EMG_prev - padding_EMG_prev;
                    y_max_dyn_EMG_prev = y_max_EMG_prev + padding_EMG_prev;
                end
                if y_min_dyn_EMG_prev >= y_max_dyn_EMG_prev
                    y_max_dyn_EMG_prev = y_min_dyn_EMG_prev + 0.1;
                end
                ylim(G.A7a, [y_min_dyn_EMG_prev, y_max_dyn_EMG_prev]);
                set(G.A7a, 'XLimMode','manual', 'YLimMode','manual');
                
                % Plot the previous EEG data (in volts)
                line(G.A7a, t_prev, curEMG_prev, 'Color','k', 'LineWidth', 1);
                set(G.A7a, 'YLim', [stableMin_EMG, stableMax_EMG], 'YLimMode', 'manual');
    
                
                % --- Dynamically Set YTick Values for g.A6a ---
                % desiredNumTicks = 4;
                % tickPrecision = 1;
                % y_min_dyn_EMG_prev_mV = y_min_dyn_EMG_prev * 1000;
                % y_max_dyn_EMG_prev_mV = y_max_dyn_EMG_prev * 1000;
                % customTicks_EMG_prev_mV = linspace(y_min_dyn_EMG_prev_mV, y_max_dyn_EMG_prev_mV, desiredNumTicks);
                % customTicks_EMG_prev_mV(1) = round(y_min_dyn_EMG_prev_mV, tickPrecision);
                % customTicks_EMG_prev_mV(end) = round(y_max_dyn_EMG_prev_mV, tickPrecision);
                % EMGprev_amp_min_mV = customTicks_EMG_prev_mV(1) * 1000;
                % EMGprev_amp_max_mV = customTicks_EMG_prev_mV(end) * 1000;
                % EMGprev_amplitudeStr = sprintf('%d to %d mV', EMGprev_amp_min_mV, EMGprev_amp_max_mV);
                % % set(G.EMGprev_ampDisplay, 'String', EMGprev_amplitudeStr);
                % if numel(unique(customTicks_EMG_prev_mV)) < desiredNumTicks || any(diff(customTicks_EMG_prev_mV) <= 0)
                %     delta = (y_max_dyn_EMG_prev_mV - y_min_dyn_EMG_prev_mV) / (desiredNumTicks - 1);
                %     customTicks_EMG_prev_mV = y_min_dyn_EMG_prev_mV : delta : y_max_dyn_EMG_prev_mV;
                %     if numel(customTicks_EMG_prev_mV) ~= desiredNumTicks
                %         customTicks_EMG_prev_mV = linspace(y_min_dyn_EMG_prev_mV, y_max_dyn_EMG_prev_mV, desiredNumTicks);
                %     end
                % end
                % customTicks_EMG_prev_mV = round(customTicks_EMG_prev_mV, tickPrecision);
                % set(G.A7a, 'YTick', customTicks_EMG_prev_mV / 1000, 'YTickMode', 'manual');
                % % For display, create formatted tick labels (in mV)
                % formattedTicks_EMG_prev = arrayfun(@(x) sprintf('%.1f', x), customTicks_EMG_prev_mV, 'UniformOutput', false);
                % set(G.A7a, 'YTickLabel', formattedTicks_EMG_prev, 'YTickMode', 'manual');
                % set(G.A7a,'XTick',[])
                % 
                % 
                % 
                % 
                % % Verification output for previous minute:
                % currentYLim_EMG_prev = get(G.A7a, 'YLim');
                % disp(['Dynamic EEG Y-axis limits for previous minute: ', num2str(currentYLim_EMG_prev)]);
                % disp(['Underlying EEG YTick values for previous minute: ', mat2str(customTicks_EMG_prev_mV)]);
                % disp(['Displayed EEG YTick labels for previous minute: ', strjoin(formattedTicks_EMG_prev, ', ')]);
                % 
    
            end
            
         
            % Plot Progress Button
            tp = G.timepointH; % time in seconds at the center of the screen
            if G.index < G.mid
                tp = gi*G.epochLen/3600-G.epochLen/3600/2;
            end
            if G.nbins - G.index < (G.mid-1)
                tp = gi*G.epochLen/3600-G.epochLen/3600/2;
            end
            li = get(G.A2,'xlim');
            cla(G.A2);
            hold(G.A2,'on')
            xlim(G.A2,li);
            set(G.A2,'YTick',[],'XTick',[],'XLimMode','manual', 'YLimMode','manual');
            
            % unless we're at the beginning or end
            if G.index < G.mid  || G.nbins - G.index < (G.mid-1)
                plot(G.A2,G.timepointH, 0.5, 'rd', 'LineWidth', 3,'MarkerFaceColor','r');
                if G.index <= (G.mid-1)
                    plot(G.A2,[0, G.epochLen/3600*G.show], [0.5,0.5], 'r','LineWidth',2);
                else
                    plot(G.A2,[G.specTh(end-G.show)+G.epochLen/3600/2, G.specTh(end)+G.epochLen/3600/2],...
                        [0.5,0.5], 'r','LineWidth',2);
                end
            else
                plot(G.A2,G.timepointH, 0.5, 'rd', 'LineWidth', 3,'MarkerFaceColor','r');
                line(G.A2,[tp-G.epochLen/3600*(G.show/2),tp+G.epochLen/3600*(G.show/2)], [0.5,0.5],...
                    'Color','r','LineWidth',2);
            end
            
            
            if tp<(li(1)+.35*diff(li)) && li(1) > G.lims(1) % we are far to the left
                xlim(G.A2,li - min([li(1)-G.lims(1), li(1)+.35*diff(li)-tp]))
            else
                if tp>(li(1)+.65*diff(li)) && li(2) < G.lims(2) % far to the right
                    xlim(G.A2,li + min([G.lims(2)-li(2), tp-li(1)-.65*diff(li)]))
                end
            end
        end
    
        function updateState() % update the sleep stage image
            
            li=xlim(G.A1);
            cla(G.A1);
            hold(G.A1, 'on');
            box(G.A1,'on');
            ylim(G.A1,[0 1]);
            ylim(G.A1,[.5 3.5]);
            xlim(G.A1,li) % make sure x limits are correct
            set(G.A1,'XLimMode','manual','YLimMode','manual');
            imagesc(G.A1,G.specTh,[1 2 3],makeSleepStageImage(G.labels),[0 4]);
            colormap(G.A1,G.colors);
            set(G.A1, 'XTickLabel', [],'XTick',[], 'YTick', [1 2 3], 'YTickLabel', {'REM', 'Wake', 'NREM'});
        end
    
        function [im] = makeSleepStageImage(state) % create the image to show
            % in the top sleep stage panel
            im = zeros(3,length(state));
            for i = 1:3
                im(i,:) = (state==i).*i;
                im(i,state==4) = 4;
            end
        end
    
    % Process keypresses
        function keypress(~, evt)
            
            switch evt.Key
                case {'rightarrow', 'uparrow'} % advance one time step
                    if G.index < G.nbins
                        G.index = G.index + 1;
                        G.timepointS  = G.specTs(G.index);
                        G.timepointH  = G.specTh(G.index);
                        updatePlots;
                    end
                    
                case {'leftarrow', 'downarrow'} % move back one time step
                    if G.index > 1
                        G.index = G.index - 1;
                        G.timepointS = G.specTs(G.index);
                        G.timepointH = G.specTh(G.index);
                        updatePlots;
                    end
                    
                case 'pageup' % jump to next bin with undefined state
                    idx = find(G.labels==4);
                    if ~isempty(idx) && any (idx > G.index)
                        idx = idx(idx>G.index);
                        G.index = idx(1);
                        G.timepointS = G.specTs(G.index);
                        G.timepointH = G.specTh(G.index);
                        updatePlots;
                    end
                    
                case 'pagedown' % jump to previous bin with undefined state
                    idx = find(G.labels==4);
                    if ~isempty(idx) && any (idx < G.index)
                        idx = idx(idx<G.index);
                        G.index = idx(end);
                        G.timepointS = G.specTs(G.index);
                        G.timepointH = G.specTh(G.index);
                        updatePlots;
                    end
                    
                case 'space' % jump to next bin with different state than current bin
                    idx = find(G.labels~=G.labels(G.index));
                    if ~isempty(idx) && any (idx > G.index)
                        idx = idx(idx>G.index);
                        G.index = idx(1);
                        G.timepointS = G.specTs(G.index);
                        G.timepointH = G.specTh(G.index);
                        updatePlots;
                    end
                    
                case 'insert' % toggle auto-scroll mode
                    G.advance = ~G.advance;
                    set(G.autobox,'Value',G.advance);
                    
                case 'a' % jump to point on spectrogram
                    axes(G.A1);
                    [x, ~] = ginput(1);
                    [G.timepointH, G.index] = findTime(G.specTh, x);
                    G.timepointS = G.specTs(G.index);
                    updatePlots;
                    
                case 'add' % zoom in
                    axes(G.A3);
                    curlims = xlim;
                    xlim([max(curlims(1), G.timepointH-.45*diff(curlims)) min(curlims(2),...
                        G.timepointH+.45*diff(curlims))]);
                    
                case 'subtract' % zoom out
                    axes(G.A3);
                    curlims = xlim;
                    xlim([max(G.lims(1), G.timepointH-1.017*diff(curlims)) min(G.lims(2),...
                        G.timepointH+1.017*diff(curlims))]);
                    
                case 'numpad0' % reset zoom level
                    axes(G.A3);
                    xlim(G.lims);
                    
                case {'r','1','numpad1'} % set to REM
                    G.labels(G.index) = 1;
                    G.unsavedChanges = 1;
                    updateState;
                    updatePlots;
                    advance;
                    
                case {'w','2','numpad2'} % set to wake
                    G.labels(G.index) = 2;
                    G.unsavedChanges = 1;
                    updateState;
                    updatePlots;
                    advance;
                    
                case {'s','3','numpad3'} % set to NREM
                    G.labels(G.index) = 3;
                    G.unsavedChanges = 1;
                    updateState;
                    updatePlots;
                    advance;
                    
                case {'x','4','numpad4'} % set to undefined
                    G.labels(G.index) = 4;
                    G.unsavedChanges = 1;
                    updateState;
                    updatePlots;
                    advance;
                    
                case 'f' % save file
                    saveFile();
                    
                case 'h' % show help menu
                    showHelp(G.A1, []);
                    
                case 'multiply' % apply label to range of timepoints
                    t = text(G.A5,18.75,-.9,sprintf(['Move the ROI\n',...
                        'boundaries,\nand then\ndouble-click it',...
                        '\nor press\nescape']),'Color','r');
                    axes(G.A1);
                    set(G.A1,'Clipping','off')
                    xl = xlim(G.A1);
                    d = diff(xl);
                    roi = imrect(G.A1,[xl(1)+d/2 - d/24,-4.0287,d/12,4.9080]);
                    rectPosition = round(wait(roi)./(G.epochLen/3600));
                    roi.delete();
                    set(G.A1,'Clipping','on')
                    if isempty(rectPosition)
                        delete(t);
                        return
                    end
                    
                    [label,~] = listdlg('PromptString','Set label to:',...
                        'SelectionMode','single',...
                        'ListString',{'REM', 'Wake','NREM','Undefined'});
                    if isempty(label) % no label selected
                        delete(t);
                        return
                    end
                    
                    idx1 = max([1,rectPosition(1)]); % starting index
                    idx2 = min([G.nbins,rectPosition(1)+rectPosition(3)]); % ending index
                    G.labels(idx1 : idx2) = label;
                    G.unsavedChanges = 1;
                    delete(t); % remove helper text
                    updatePlots; % update plots
                    updateState;
            end
        end
    
    % functions called by button presses
        function brightSpect(src,~)
            
            G.cmax = G.cmax - G.cmax/10;
            G.caxis1 = [G.caxis1(1), G.cmax];
            caxis(G.A3,G.caxis1)
            updatePlots;
            defocus(src);
        end
    
        function dimSpect(src,~)
            
            G.cmax = G.cmax + G.cmax/10;
            G.caxis1 = [G.caxis1(1), G.cmax];
            caxis(G.A3,G.caxis1)
            updatePlots;
            defocus(src);
        end
    
        function setState(src,~, a) % apply label to single time bin
            evt= struct;
            keys = {'r','w','s'};
            evt.Key = keys{a};
            keypress([], evt);
            defocus(src);
        end
    
        function setRange(src, ~) % apply label to range of time bins
            evt= struct;
            evt.Key = 'multiply';
            keypress([], evt);
            defocus(src);
        end
    
        function fct_showmenu(src, ~)
            
           options = [1, 3, 5, 7, 9, 15];
            G.show = options(src.Value);
            G.mid = ceil(G.show/2);
            updatePlots;
            defocus(src);
        end
    
        function fct_shiftupEEG(src,~)
            
            G.eegYlim = [G.eegYlim(1)+.04*diff(G.eegYlim), G.eegYlim(2)+.04*diff(G.eegYlim)];
            updatePlots;
            defocus(src);
        end
    
        function fct_shiftdownEEG(src,~)
            
            G.eegYlim = [G.eegYlim(1)-.04*diff(G.eegYlim), G.eegYlim(2)-.04*diff(G.eegYlim)];
            updatePlots;
            defocus(src);
        end
    
        function fct_shiftupEMG(src,~)
            
            G.emgYlim = [G.emgYlim(1)+.04*diff(G.emgYlim), G.emgYlim(2)+.04*diff(G.emgYlim)];
            updatePlots;
            defocus(src);
        end
    
        function fct_shiftdownEMG(src,~)
            
            G.emgYlim = [G.emgYlim(1)-.04*diff(G.emgYlim), G.emgYlim(2)-.04*diff(G.emgYlim)];
            updatePlots;
            defocus(src);
        end
    
        function fct_selectlocation(src,~)
            evt= struct;
            evt.Key = 'a';
            keypress([], evt);
            defocus(src);
        end
    
        function fct_zoomin_t(src,~)
            
            axes(G.A3);
            curlims = xlim;
            xlim([max(curlims(1), G.timepointH-.45*diff(curlims)) min(curlims(2), G.timepointH+.45*diff(curlims))]);
            defocus(src);
        end
    
        function fct_zoomout_t(src,~)
            
            axes(G.A3);
            curlims = xlim;
            xlim([max(G.lims(1), G.timepointH-1.02*diff(curlims)) min(G.lims(2), G.timepointH+1.02*diff(curlims))]);
            defocus(src);
        end
    
        function fct_zoominEEG(src, ~)
            % Zoom in on the EEG axis (G.A7) by reducing its y-range by 5% on each side.
            % Print statements are added to check the current and new YLim values.
            
            % Retrieve current y-axis limits (in volts)
            current_EEGYLim = get(G.A7, 'YLim');
            fprintf('Before zoom: current YLim = [%.6f, %.6f] V\n', current_EEGYLim(1), current_EEGYLim(2));
            
            % Calculate the current range
            rangeY = diff(current_EEGYLim);
            fprintf('Current range = %.6f V\n', rangeY);
            
            % Compute new y-axis limits by reducing the range by 5% on both sides
            newYLim = [current_EEGYLim(1) + 0.15 * rangeY, current_EEGYLim(2) - 0.15 * rangeY];
            fprintf('New YLim = [%.6f, %.6f] V\n', newYLim(1), newYLim(2));
            
            % Apply the new y-axis limits, and keep them fixed
            set(G.A7, 'YLim', newYLim, 'YLimMode', 'manual');
            
            % Update the plots and remove focus from the button
            updatePlots;
            defocus(src);
        end
    
        function fct_zoomoutEEG(src,~)
            
            G.eegYlim = [G.eegYlim(1)-.05*diff(G.eegYlim), G.eegYlim(2)+.05*diff(G.eegYlim)];
            updatePlots;
            defocus(src);
        end
    
        function fct_zoominEMG(src,~)
            
            G.emgYlim = [G.emgYlim(1)+.05*diff(G.emgYlim), G.emgYlim(2)-.05*diff(G.emgYlim)];
            updatePlots;
            defocus(src);
        end
    
        function fct_zoomoutEMG(src,~)
            
            G.emgYlim = [G.emgYlim(1)-.05*diff(G.emgYlim), G.emgYlim(2)+.05*diff(G.emgYlim)];
            updatePlots;
            defocus(src);
        end
    
        function fct_zoomreset_t(src,~) % reset zoom level
            
            axes(G.A3);
            xlim(G.lims);
            defocus(src);
        end
    
        function scrollCallback(a,~) % respond to user input in the auto-scroll box
            
            G.advance = a.Value;
            defocus(a);
        end
    
        function showHelp(src,~) % show help menu
            createmode = struct;
            createmode.WindowStyle = 'non-modal';
            createmode.Interpreter = 'tex';
            msgbox({'\fontsize{13}\fontname{Courier New}User manual:';
                'After closing this window, click inside the figure to resume.';
                'The lower panels are a zoomed-in subset of data in the upper panels.';
                'The red diamond on top of the EEG shows the current time point,';
                'and the red lines extending on either side mark the subset shown ';
                'in the lower panels.';
                ' ';
                'Keyboard shortcuts:';
                'Up/right arrow : scroll one time step forward';
                'Down/left arrow : scroll one time step backward';
                'A : choose time point: click to jump to a location';
                '+/- : zoom in/out';
                '0 (zero) : reset zoom';
                'F : save sleep stage labels - see code for details';
                'pg up/down : skip to next/previous undefined period';
                'space : skip to first bin with different state than current bin';
                'insert : toggle auto-scroll mode, which advances one time step';
                '         automatically after applying a label';
                'H : show help menu';
                ' ';
                'R or 1 : set label as REM (cyan)';
                'W or 2 : Wake (magenta)';
                'S or 3 : NREM (gray)';
                'X or 4 : Undefined (black)';
                '* (multiply) : set labels for a range of timepoints.';
                'On the upper sleep stage panel, adjust the ROI to';
                'select the timepoints, then double-click it.';
                'Then, choose the sleep stage from the menu.';
                }, createmode);
            defocus(src);
        end
    
        function saveCallback(src,~)
            saveFile();
            defocus(src);
        end
    
        function loadFile(src,~) % load sleep stage labels from file
            
            [file,path] = uigetfile('*.mat');
            if ~ischar(file)
                msgbox('No file specified, file not loaded')
                return
            end
            f = load([path,file]); % load the file
            if checkLabels(f,1) % make sure file has the correct contents
                G.labels = f.labels; % use these labels
                G.savepath = ''; % clear the save path to avoid confusion
                updateState; % plot the new labels
            end
            defocus(src);
        end
    
        function defocus(src) % remove focus from a clicked button so that keypresses work again
            set(src, 'Enable', 'off');
            drawnow;
            set(src, 'Enable', 'on');
        end
    
    % other functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function saveFile()
            
            if isempty(G.savepath) % get file path if we need it
                [file,path] = uiputfile('*.mat');
                if ~ischar(file) % if no file given
                    msgbox('No file specified, file not saved')
                    return
                end
                G.savepath = [path,file]; % store this filename
            end
            labels = G.labels; % store sleep stage labels in a variable called 'labels'
            save(G.savepath, 'labels'); % write labels to file
            G.unsavedChanges = 0;
        end
    
    % check if a structure x has a properly constructed labels field
        function [allChecksPassed] = checkLabels(x, checkLen)
            allChecksPassed = 0;
            if ~isfield(x, 'labels') % needs a labels variable
                msgbox('Error: file does not contain a variable called labels')
                return
            end
            if ~isa(x.labels,'numeric') % labels are numeric
                msgbox('Error: labels variable must be numeric')
                return
            end
            if ~isempty(setdiff(unique(x.labels), 1:4)) % in the range 1:4
                msgbox('Error: labels variable must be in the range 1:4')
                return
            end
            if length(size(x.labels)) > 2 || min(size(x.labels)) >  1 % should be a vector
                msgbox('Error: labels variable must be a vector')
                return
            end
            if checkLen
                
                if length(x.labels) ~= G.nbins % in the range 1:4
                    msgbox(['Error: labels must be of length ',num2str(G.nbins)])
                    return
                end
            end
            if isrow(x.labels) % ensure that it's a column
                x.labels = x.labels';
            end
            allChecksPassed = 1;
        end
    
        function [timeString] = sec2hr(s) % convert seconds into hr:mn:sec.abcd
            m = floor(s/60);
            s2 = mod(s,60);
            h = floor(m/60);
            m2 = mod(m,60);
            timeString = sprintf('%02d:%02d:%05.2f',h,m2,s2);
        end
    
        function [t, index] = findTime(time, x) % find closest value of time to x
            [~, index] = min(abs(time - x));
            t = time(index);
        end
    
        function advance() % if auto-scrolling is on, advance to the next time bin
            if G.advance
                if G.index < G.nbins
                    G.index = G.index + 1;
                    G.timepointS  = G.specTs(G.index);
                    G.timepointH  = G.specTh(G.index);
                    if G.index < G.nbins-1 && G.index > 2
                        updatePlots;
                    end
                end
            end
        end
    
        function closeReq(~,~) %
            if G.unsavedChanges
                answer = questdlg('There are unsaved changes. Really quit?', ...
                    'Unsaved changes', ...
                    'Quit','Cancel','Cancel');
                if ~strcmp(answer,'Quit')
                    return
                end
            end
            delete(gcf)
        end
    
    end
    