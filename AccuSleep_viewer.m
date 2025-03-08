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
G =     struct; % holds everything

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
emgCap = mean(G.cappedEMG) + 3*std(G.cappedEMG);
G.cappedEMG(G.cappedEMG > emgCap) = emgCap;

% set various parameters
G.show = 15; %  default number of bins to display on screen
G.dt = 1/G.originalSR; % duration of each EEG/EMG sample in seconds
G.advance = 1; % whether to advance automatically when a state is assigned
G.colors = [1 1 1; .47 .67 .19; .14 .2 .57; 0.996 0.758 0.039; 0.25 0.45 0.05; 0.05 0.15 0.50;0.85 0.65 0.03]; %colors for sleep stages
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
else % if no sleep stages provided, set to W*
    G.labels = ones(G.nbins,1) * 5;
end

% get spectrogram and time axes
showFreqs = find(fAxis <= 21); % only show frequencies under 25 Hz
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
G.A5 = axes('Units', 'Normalized', 'Position', [0 .7 .05 .02],'XColor','w','YColor','w');
hold(G.A5,'on');
set(G.A5,'Box','off','XLim',[0 1],'YLim',[0 0.1],'XTick',[],'YTick',[],'Clipping','off')
G.A5.Toolbar.Visible = 'off';

% EEG
% spectrogram
G.A3 = axes('Units', 'Normalized', 'Position', [0.05 0.935 0.87 0.06]);
G.A3.Toolbar.Visible = 'off';
set(gca, 'FontSize', 10, 'LineWidth', 2, 'XTick', [], 'YTick', []);
ylabel(G.A3, 'Spec.');

% Upper time point indicator
G.A8 = axes('Units', 'Normalized', 'Position', [0.05 0.92  0.87 0.015],'XTick',[],'YTick',[]);
G.A8.Toolbar.Visible = 'off';

% Upper sleep stage labels
G.A1 = axes('Units', 'Normalized', 'Position', [0.05 0.865 0.87 0.05]);
G.A1.Toolbar.Visible = 'off';

% EEG signal - 1 Minute
G.A6a = axes('Units', 'Normalized', 'Position', [0.05 0.78, 0.87 .09]);
G.A6a.Toolbar.Visible = 'off';
ylabel(G.A6a, 'EEG-1');
set(G.A6a,'XTick',[])

% % EMG signal - 1 Minute
G.A7a = axes('Units', 'Normalized', 'Position', [0.05 0.68 0.87 .09]);
G.A7a.Toolbar.Visible = 'off';
ylabel(G.A7a, 'EMG-1');
set(G.A7a,'XTick',[])

% EEG signal  - main
G.A6 = axes('Units', 'Normalized', 'Position', [0.05 0.585 0.87 .09]);
G.A6.Toolbar.Visible = 'off';
ylabel(G.A6, 'EEG');

% EMG signal  - main
G.A7 = axes('Units', 'Normalized', 'Position', [0.05 0.475 0.87 .09]);
G.A7.Toolbar.Visible = 'off';
ylabel(G.A7, 'EMG');
set(G.A7,'XTick',[])

% EEG signal + 1 minute
G.A6b = axes('Units', 'Normalized', 'Position', [0.05 0.375 0.87 .09]);
G.A6b.Toolbar.Visible = 'off';
ylabel(G.A6b, 'EEG+1');
set(G.A6b,'XTick',[])

% EMG signal + 1 minute
G.A7b = axes('Units', 'Normalized', 'Position', [0.05 0.275 0.87 .09]);
G.A7b.Toolbar.Visible = 'off';
ylabel(G.A7b, 'EMG+1');
set(G.A7b,'XTick',[])

% Current PSD
G.A10 = axes('Units', 'Normalized', 'Position', [0.44 0.155 0.085 0.10]);
G.A10.Toolbar.Visible = 'off';


% -4 PSD
G.A17 = axes('Units', 'Normalized', 'Position', [0.05 0.155 0.085 0.10]);
G.A17.Toolbar.Visible = 'off';
% -3 PSD
G.A18 = axes('Units', 'Normalized', 'Position', [0.145 0.155 0.085 0.10]);
G.A18.Toolbar.Visible = 'off';
% -2 PSD
G.A19 = axes('Units', 'Normalized', 'Position', [0.24 0.155 0.085 0.10]);
G.A19.Toolbar.Visible = 'off';
% -1 PSD
G.A20 = axes('Units', 'Normalized', 'Position', [0.34 0.155 0.085 0.10]);
G.A20.Toolbar.Visible = 'off';
% +1 PSD
G.A21 = axes('Units', 'Normalized', 'Position', [0.54 0.155 0.085 0.10]);
G.A21.Toolbar.Visible = 'off';
% +2 PSD
G.A22 = axes('Units', 'Normalized', 'Position', [0.645 0.155 0.085 0.10]);
G.A22.Toolbar.Visible = 'off';
% +3 PSD
G.A23 = axes('Units', 'Normalized', 'Position', [0.74 0.155 0.085 0.10]);
G.A23.Toolbar.Visible = 'off';
% +4 PSD
G.A24 = axes('Units', 'Normalized', 'Position', [0.84 0.155 0.085 0.10]);
G.A24.Toolbar.Visible = 'off';


% Lower sleep stage labels
G.A9 = axes('Units', 'Normalized', 'Position', [0.05 0.01 0.87 0.1]);

% Plot the EMG data (in volts)
globalMin_EMG = min(G.EMG);
globalMax_EMG = max(G.EMG);
globalRange_EMG = globalMax_EMG - globalMin_EMG;
padding_EMG = 0.02 * globalRange_EMG;
% fprintf('Global EMG Range = %.8f , %.8f \n', globalMax_EMG, globalMin_EMG);
stableMin_EMG = globalMin_EMG - padding_EMG;
stableMax_EMG = globalMax_EMG + padding_EMG;
% --- In updatePlots, apply the stable EEG axis limits ---
set(G.A7, 'YLim', [stableMin_EMG, stableMax_EMG], 'YLimMode', 'manual');
set(G.A7a, 'YLim', [stableMin_EMG, stableMax_EMG], 'YLimMode', 'manual');
set(G.A7b, 'YLim', [stableMin_EMG, stableMax_EMG], 'YLimMode', 'manual');




% --- Compute Global EMG Limits (in volts) ---
globalMin_EEG = min(G.EEG);
globalMax_EEG = max(G.EEG);
globalRange_EEG = globalMax_EEG - globalMin_EEG;
padding_EEG = 0.02 * globalRange_EEG;
% fprintf('Global EEG Range = %.8f , %.8f \n', globalMax_EEG, globalMin_EEG);

stableMin_EEG = globalMin_EEG - padding_EEG;
stableMax_EEG = globalMax_EEG + padding_EEG;

% --- In updatePlots, apply the stable EEG axis limits ---
set(G.A6, 'YLim', [stableMin_EEG, stableMax_EEG], 'YLimMode', 'manual');
set(G.A6a, 'YLim', [stableMin_EEG, stableMax_EEG], 'YLimMode', 'manual');
set(G.A6b, 'YLim', [stableMin_EEG, stableMax_EEG], 'YLimMode', 'manual');

linkaxes([G.A1,G.A3, G.A8], 'x'); % upper panel x axes should stay linked

% buttons
G.helpbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized','BackgroundColor',[1 .8 .8],...
    'Position',[.93 .01 .062 .055],'String','Help','Callback',@showHelp,'FontSize',9,...
    'ToolTip','Show help menu (H)');
G.savebtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized','BackgroundColor',[.8 1 .8],...
    'Position',[.93 .07 .062 .045],'String','Save labels','Callback',@saveCallback,'FontSize',9,...
    'ToolTip','Save labels to file (F)');
G.loadbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
    'Position',[.93 .12 .062 .03],'String','Load labels','Callback',@loadFile,'FontSize',9,...
    'ToolTip','Load labels from file');
G.brightbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
    'Position',[.93 .88 .062 .025],'String','Brighter','Callback',@brightSpect,...
    'FontSize',9,'ToolTip','Make EEG spectrogram brighter');
G.dimbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
    'Position',[.93 .85 .062 .025],'String','Dimmer','Callback',@dimSpect,...
    'FontSize',9,'ToolTip','Make EEG spectrogram dimmer');
G.selectbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
    'Position',[.93 .78 .062 .04],'String','<html>Select<br>timepoint',...
    'Callback',@fct_selectlocation,'FontSize',9,'ToolTip','Click a timepoint to show (A)');
G.zoominbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
    'Position',[.93 .97 .062 .025],'String','Zoom IN','Callback',@fct_zoomin_t,...
    'FontSize',9,'ToolTip','Increase zoom level (+)');
G.zoomoutbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
    'Position',[.93 .94 .062 .025],'String','Zoom OUT','Callback',@fct_zoomout_t,...
    'FontSize',9,'ToolTip','Decrease zoom level (-)');
G.zoomresetbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
    'Position',[.93 .91 .062 .025],'String','Reset zoom','Callback',@fct_zoomreset_t,...
    'FontSize',9,'ToolTip','Reset zoom level (0)');
G.gui_zoominEEG = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
    'Position',[0.93 0.555 0.03 0.03],'Callback', @fct_zoominEEG,'String','EEG+',...
    'FontSize',9);
G.gui_zoomoutEEG = uicontrol(WIN,'Style','pushbutton','Units','normalized',...
    'Position',[0.93 0.515 0.03 0.03],'Callback', @fct_zoomoutEEG,'String','EEG-',...
    'FontSize',9);
G.gui_zoominEMG = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
    'Position',[0.93 0.41 0.03 0.03],'Callback', @fct_zoominEMG,'String','EMG+',...
    'FontSize',9);
G.gui_zoomoutEMG = uicontrol(WIN,'Style','pushbutton','Units','normalized',...
    'Position',[0.93 0.37 0.03 0.03],'Callback', @fct_zoomoutEMG,'String','EMG-',...
    'FontSize',9);
G.showMenu = uicontrol(WIN,'Style','popupmenu','Units','normalized',...
    'Position',[0.93 0.75 0.062 0.02],'Callback', @fct_showmenu,...
    'String',{'1 epoch','3 epochs','5 epochs','7 epochs','9 epochs','15 epochs'},...
    'Value',6);

G.meanDelta = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [.925 .24 .07 .035], 'String', '', ...
    'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
G.meanTheta = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [.925 .20 .07 .035], 'String', '', ...
    'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
G.thetaRatio = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [.925 .16 .07 .035], 'String', '', ...
    'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
G.eegIntegral = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [.925 .32 .07 .035], 'String', '', ...
    'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
G.emgIntegral = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [.925 .28 .07 .035], 'String', '', ...
    'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
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
set(G.A3, 'YTick', 0:5:20, 'YTickLabel', 0:5:20);
ylabel(G.A3, 'Spec.');


updateState;

% EMG_integral = computeSignalIntegral(G.processedEMG,G.SR);
% G.EMG_Integral = EMG_integral;

[totalfreq, totalpsdVals] = computeBinPSD(G.EEG, G.SR, 4, 0.9);
meanVal = mean(totalpsdVals);
stdDev = std(totalpsdVals);

ylim(G.A10, [0, meanVal + 3*stdDev]);
ylim(G.A20, [0, meanVal + 3*stdDev]);
ylim(G.A21, [0, meanVal + 3*stdDev]);


% fprintf("Total PSD Vals: %s \n %s",min(totalpsdVals), max(totalpsdVals))
% Plot everything else
updatePlots;


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

        seq=G.labels((1:G.show)+gi-G.mid+(G.mid-gi)*(gi<G.mid)-...
            (gi-G.nbins+(G.mid-1))*(gi>(G.nbins-(G.mid-1))));
        x = -n:n;
        cla(G.A9)
        hold(G.A9,'on')
        xlim(G.A9, [-n-0.5 n+0.5]);
        ylim(G.A9, [0.5 3.5]);
        set(G.A9, 'XLimMode','manual', 'YLimMode','manual');
        for i = 1:length(seq)
            if seq(i) == 4 || seq(i) == 5 || seq(i) == 6
                % For stages 4, 5, and 6, draw a large block
                pX = [x(i)-0.5, x(i)+0.5, x(i)+0.5, x(i)-0.5];
                pY = [3.5, 3.5, 0.5, 0.5];
                patch(G.A9, pX, pY, G.colors(seq(i)+1, :), 'EdgeColor', 'none');
            else
                % For other stages, draw a smaller block
                pY = [seq(i)+0.5, seq(i)+0.5, seq(i)-0.5, seq(i)-0.5];
                pX = [x(i)-0.5, x(i)+0.5, x(i)+0.5, x(i)-0.5];
                patch(G.A9, pX, pY, G.colors(seq(i)+1, :), 'EdgeColor', 'none');
            end
        end

        set(G.A9, 'XTickLabel', [],'XTick',[], 'YTick', [1 2 3], 'YTickLabel', {'REM', 'Wake', 'NREM'});

        line(G.A9, [-0.5, -0.5], [0.5,1], 'Color','r','LineWidth',1);
        line(G.A9, [0.5, 0.5], [0.5,1], 'Color','r','LineWidth',1);

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
        % curEMGIntegral = G.EMG_Integral(ii);
        y_min_EMG = min(curEMG);
        y_max_EMG = max(curEMG);
        yr_EMG = y_max_EMG - y_min_EMG;
        padding_EMG = 0.02 * yr_EMG;

        % set emgIntegral
        % set(G.emgIntegral, 'String', sprintf('EMG Integral: %.10f ', curEMGIntegral));


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

        line(G.A7, t, curEMG, 'Color', 'k', 'LineWidth', 1);

        currentYLim = get(G.A7, 'YLim');  % currentYLim = [bottom, top]
        % Define an indicator height as a fraction of the current range
        indicatorHeight_EMG = 0.125 * diff(currentYLim);

        % Draw vertical lines at the left and right epoch boundaries, starting at the bottom (currentYLim(1))
        line(G.A7, ones(1,2)*(G.timepointS - G.epochLen/2), [currentYLim(1), currentYLim(1) + indicatorHeight_EMG], 'Color', 'r', 'LineWidth', 0.5);
        line(G.A7, ones(1,2)*(G.timepointS + G.epochLen/2), [currentYLim(1), currentYLim(1) + indicatorHeight_EMG], 'Color', 'r', 'LineWidth', 0.5);

        % Draw a horizontal line connecting the tops of the vertical markers at the bottom
        line(G.A7, [G.timepointS - G.epochLen/2, G.timepointS + G.epochLen/2], [currentYLim(1), currentYLim(1)], 'Color', 'r', 'LineWidth', 0.5);

        cla(G.A6)
        hold(G.A6, 'on');
        xlim(G.A6,[t(1)-G.dt t(end)]);
        % Compute dynamic y-axis limits based on current EEG window data (in volts)

        curEEG = G.EEG(ii);
        y_min_EEG = min(curEEG);
        y_max_EEG = max(curEEG);
        yr_EEG = y_max_EEG - y_min_EEG;
        padding_EEG = 0.02 * yr_EEG;
        % Compute the integral of the current EMG
        % curEEG_integral = computeSignalIntegral(curEEG,G.SR);
        % set(G.eegIntegral, 'String', sprintf('EEG Integral: %.10f ', curEEG_integral));
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

        % Draw vertical red lines at the left and right boundaries of the current epoch,
        % with the line extending from the top (y_max_dyn_EEG) down by indicatorHeight_EEG.
        currentYLim = get(G.A6, 'YLim');  % currentYLim = [bottom, top]
        % Define an indicator height as a fraction of the current range
        indicatorHeight_EMG = 0.125 * diff(currentYLim);

        % Draw vertical lines at the left and right epoch boundaries, starting at the bottom (currentYLim(1))
        line(G.A6, ones(1,2)*(G.timepointS - G.epochLen/2), [currentYLim(2), currentYLim(2) - indicatorHeight_EMG], 'Color', 'r', 'LineWidth', 0.5);
        line(G.A6, ones(1,2)*(G.timepointS + G.epochLen/2), [currentYLim(2), currentYLim(2) - indicatorHeight_EMG], 'Color', 'r', 'LineWidth', 0.5);

        % Draw a horizontal line connecting the tops of the vertical markers at the bottom
        line(G.A6, [G.timepointS - G.epochLen/2, G.timepointS + G.epochLen/2], [currentYLim(2), currentYLim(2)], 'Color', 'r', 'LineWidth', 0.5);

        % label x axis nicely
        G.A6.XTick = tp-(G.show/2)*G.epochLen + G.epochLen*(0:G.show);
        ticks = G.A6.XTick;
        xlbl = cell(1, length(ticks));
        for i = 1:length(ticks)
            xlbl{i} = sec2hr(ticks(i));
        end
        G.A6.XTickLabel = xlbl;

        % Plot Progress Button
        tp = G.timepointH; % time in seconds at the center of the screen
        if G.index < G.mid
            tp = gi*G.epochLen/3600-G.epochLen/3600/2;
        end
        if G.nbins - G.index < (G.mid-1)
            tp = gi*G.epochLen/3600-G.epochLen/3600/2;
        end

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


        showAllEpochsPSD(G)

                %% Plot EEG signal - Previous Minute in zoomineeg
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

            % Plot the previous EEG data (in volts)
            line(G.A6a, t_prev, curEEG_prev, 'Color','k', 'LineWidth', 1);

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
            % ylim(G.A7a, [y_min_dyn_EMG_prev, y_max_dyn_EMG_prev]);
            set(G.A7a, 'XLimMode','manual', 'YLimMode','manual');

            % Plot the previous EEG data (in volts)
            line(G.A7a, t_prev, curEMG_prev, 'Color','k', 'LineWidth', 1);


        end

        if (G.timepointS + 60) * G.dt <= G.eegLen * G.dt  % or simply: if G.timepointS+60 <= totalDuration
            % Define center time for the future minute:
            tp_future = G.timepointS + 60;

            % Determine window size for the future minute:
            n_future = round((G.show * G.epochLen) / G.dt / 2);
            i_center_future = round(tp_future / G.dt);
            ii_future = i_center_future - n_future : i_center_future + n_future;
            ii_future(ii_future <= 0) = 1;
            ii_future(ii_future >= G.eegLen) = G.eegLen;

            % Create time vector for the future window:
            t_future = (tp_future - n_future * G.dt) : G.dt : (tp_future + n_future * G.dt);

            % Clear and prepare the future EEG axis (G.A6b)
            cla(G.A6b);
            hold(G.A6b, 'on');
            xlim(G.A6b, [t_future(1)-G.dt, t_future(end)]);

            % Compute dynamic y-axis limits for the future EEG window (in volts)
            curEEG_future = G.EEG(ii_future);
            y_min_EEG_future = min(curEEG_future);
            y_max_EEG_future = max(curEEG_future);
            yr_EEG_future = y_max_EEG_future - y_min_EEG_future;
            padding_EEG_future = 0.02 * yr_EEG_future;

            if yr_EEG_future < eps
                y_min_dyn_EEG_future = y_min_EEG_future - 0.1;
                y_max_dyn_EEG_future = y_max_EEG_future + 0.1;
            else
                y_min_dyn_EEG_future = y_min_EEG_future - padding_EEG_future;
                y_max_dyn_EEG_future = y_max_EEG_future + padding_EEG_future;
            end
            if y_min_dyn_EEG_future >= y_max_dyn_EEG_future
                y_max_dyn_EEG_future = y_min_dyn_EEG_future + 0.1;
            end


            % Plot the future EEG data (in volts)
            line(G.A6b, t_future, curEEG_future, 'Color','k', 'LineWidth', 1);


            cla(G.A7b);
            hold(G.A7b, 'on');
            xlim(G.A7b, [t_future(1)-G.dt, t_future(end)]);

            % Compute dynamic y-axis limits for the future EEG window (in volts)
            curEMG_future = G.EMG(ii_future);
            y_min_EMG_future = min(curEMG_future);
            y_max_EMG_future = max(curEMG_future);
            yr_EMG_future = y_max_EMG_future - y_min_EMG_future;
            padding_EMG_future = 0.02 * yr_EMG_future;

            if yr_EMG_future < eps
                y_min_dyn_EMG_future = y_min_EMG_future - 0.1;
                y_max_dyn_EMG_future = y_max_EMG_future + 0.1;
            else
                y_min_dyn_EMG_future = y_min_EMG_future - padding_EMG_future;
                y_max_dyn_EMG_future = y_max_EMG_future + padding_EMG_future;
            end
            if y_min_dyn_EMG_future >= y_max_dyn_EMG_future
                y_max_dyn_EMG_future = y_min_dyn_EMG_future + 0.1;
            end


            % Plot the future EEG data (in volts)
            line(G.A7b, t_future, curEMG_future, 'Color','k', 'LineWidth', 1);

        else
            cla(G.A6b);
            text(0.5, 0.5, 'No future minute available', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10);
            cla(G.A7b);
            text(0.5, 0.5, 'No future minute available', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10);
        end


        % Plot Progress Button
        tp = G.timepointH; % time in seconds at the center of the screen
        if G.index < G.mid
            tp = gi*G.epochLen/3600-G.epochLen/3600/2;
        end
        if G.nbins - G.index < (G.mid-1)
            tp = gi*G.epochLen/3600-G.epochLen/3600/2;
        end


        li = get(G.A8,'xlim');
        cla(G.A8);
        hold(G.A8,'on')
        xlim(G.A8,li);
        set(G.A8,'YTick',[],'XTick',[],'XLimMode','manual', 'YLimMode','manual');

        % unless we're at the beginning or end
        if G.index < G.mid  || G.nbins - G.index < (G.mid-1)
            plot(G.A8,G.timepointH, 0.5, 'rd', 'LineWidth', 1.5,'MarkerFaceColor','r');

            if G.index <= (G.mid-1)
                plot(G.A8,[0, G.epochLen/3600*G.show], [0.5,0.5], 'r','LineWidth',2);
            else
                plot(G.A8,[G.specTh(end-G.show)+G.epochLen/3600/2, G.specTh(end)+G.epochLen/3600/2],...
                    [0.5,0.5], 'r','LineWidth',2);

            end
        else
            plot(G.A8,G.timepointH, 0.5, 'rd', 'LineWidth', 3,'MarkerFaceColor','r');
            line(G.A8,[tp-G.epochLen/3600*(G.show/2),tp+G.epochLen/3600*(G.show/2)], [0.5,0.5],...
                'Color','r','LineWidth',2);
        end


        if tp<(li(1)+.45*diff(li)) && li(1) > G.lims(1) % we are far to the left
            xlim(G.A8,li - min([li(1)-G.lims(1), li(1)+.45*diff(li)-tp]))
        else
            if tp>(li(1)+.55*diff(li)) && li(2) < G.lims(2) % far to the right
                xlim(G.A8,li + min([G.lims(2)-li(2), tp-li(1)-.55*diff(li)]))
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
        imagesc(G.A1,G.specTh,[1 2 3],makeSleepStageImage(G.labels),[0 6]);
        colormap(G.A1,G.colors);
        set(G.A1, 'XTickLabel', [],'XTick',[], 'YTick', [1 2 3], 'YTickLabel', {'Rem', 'Wake', 'NREM'});
    end
    function integralSignal = computeSignalIntegral(signal, fs, options)
        % Default options
        if nargin < 3 || isempty(options)
            options.baselineCorrect = false;
            options.normalize = false;
            options.method = 'trapz';
        end

        % Ensure signal is a column vector
        if isrow(signal)
            signal = signal';
        end

        % Get number of samples
        nSamples = length(signal);

        % Create time vector (in seconds)
        t = (0:nSamples-1) / fs;

        % Baseline correction (remove mean if requested)
        if options.baselineCorrect
            signal = signal - mean(signal);
            fprintf('Signal baseline corrected (mean removed).\n');
        end

        % Compute integral based on method
        switch lower(options.method)
            case 'trapz'
                % Use trapezoidal numerical integration (MATLAB's trapz)
                integralSignal = trapz(t, signal);
                % Reshape to match signal length for cumulative effect
                integralSignal = cumsum(signal) * (t(2) - t(1)); % Scale by time step
            case 'cumsum'
                % Use cumulative sum (simpler but less accurate for non-uniform data)
                dt = 1 / fs;  % Time step
                integralSignal = cumsum(signal) * dt;
            otherwise
                error('Unsupported integration method. Use ''trapz'' or ''cumsum''.');
        end

        % Normalize by total time if requested
        if options.normalize
            totalTime = t(end) - t(1);
            integralSignal = integralSignal / totalTime;
            fprintf('Integral normalized by total time (%f seconds).\n', totalTime);
        end

        % Ensure output is a column vector
        if isrow(integralSignal)
            integralSignal = integralSignal';
        end

        fprintf('Integral computed with %d samples, fs = %d Hz.\n', nSamples, fs);
    end

    function printBandMeans(freq, psd_dB)
        % Define frequency bands and a corresponding label (for printing)
        bands = {
            [0, 4],      'Red';    % 0-4 Hz => red
            [5, 9],      'Blue';   % 5-9 Hz => blue
            [10, 15],    'Green';  % 10-15 Hz => green
            [16, 20],    'Gray'    % 16-20 Hz => gray
        };

        nBands = size(bands, 1);
        bandMeans = zeros(1, nBands);

        % Calculate the mean PSD (in dB) for each frequency band
        for b = 1:nBands
            range = bands{b, 1};
            idx = (freq >= range(1)) & (freq <= range(2));
            if any(idx)
                bandMeans(b) = mean(psd_dB(idx));
            else
                bandMeans(b) = NaN;
            end
        end

        % Print the mean power for each band
        % for b = 1:nBands
        %     fprintf('Mean power in %d-%d Hz (%s): %.10f dB\n', ...
        %         bands{b, 1}(1), bands{b, 1}(2), bands{b, 2}, bandMeans(b));
        % end

        % Calculate the mean power in the 0-4 Hz band
        meanPower_0_4Hz = sum(psd_dB(freq >= 0 & freq <= 4));
        meanPower_5_9Hz = sum(psd_dB(freq >=5 & freq <=9));
        tRatio = (meanPower_5_9Hz-meanPower_0_4Hz)/meanPower_0_4Hz;

        % Update the text box
        set(G.meanDelta, 'String', sprintf('Delta: %.10f ', meanPower_0_4Hz));
        set(G.meanTheta, 'String', sprintf('Theta: %.10f ', meanPower_5_9Hz));
        set(G.thetaRatio, 'String', sprintf('Theta Ratio: %.3f ', tRatio));



        % Print the mean power in the 0-4 Hz band
        % fprintf('Mean power in 0-4 Hz: %.10f dB\n', meanPower_0_4Hz);
    end
    function showAllEpochsPSD(G)
        % 1) Clear G.A10 and hold on
        cla(G.A10);
        hold(G.A10, 'on');

        % Welch analysis parameters
        windowSeconds  = 4;      % sub-window length for Welch
        overlapPercent = 0.9;    % 90% overlap
        freqRange      = [0 20]; % show 0..20 Hz
        fs             = 1 / G.dt;  % sampling rate

        % Bins to display
        n         = (G.show -1)/2;
        % PSD Bins
        centerBin = G.index;
        centerBin_MinusFour = centerBin-4;
        centerBin_MinusThree = centerBin-3;
        centerBin_MinusTwo = centerBin-2;
        centerBin_MinusOne = centerBin-1;
        centerBin_PlusOne = centerBin+1;
        centerBin_PlusTwo = centerBin+2;
        centerBin_PlusThree = centerBin+3;
        centerBin_PlusFour = centerBin+4;
        binStart  = centerBin - n;
        binEnd    = centerBin + n;
        fprintf("Bin start is: %d and Bin End is: %d\n", binStart, binEnd);

        % We'll track min/max across all bins to standardize the y-scale
        current_data = getBinData(G, centerBin);
        [freq, psdVals] = computeBinPSD(current_data, fs, windowSeconds, overlapPercent);
        localMax = max(psdVals);
        localMin = min(psdVals);

        mask = (freq>= freqRange(1)) & (freq <= freqRange(2));
        freq = freq(mask);
        psdVals = psdVals(mask);
        psd_dB = psdVals;
        printBandMeans(freq, psd_dB);
        psd_dB = psdVals;
        % 6) Color-code each frequency point in freq by 4 bands:
        %    0-4 Hz => red, 5-9 Hz => blue, 10-15 Hz => green, 16-20 Hz => gray
        cMap = zeros(length(freq), 3);  % each row = RGB
        for i = 1:length(freq)
            f = freq(i);
            if f >= 0  && f <= 4
                cMap(i,:) = [1, 0, 0];      % red
            elseif f > 4  && f <= 9
                cMap(i,:) = [0, 0, 1];      % blue
            elseif f > 9 && f <= 15
                cMap(i,:) = [0, 1, 0];      % green
            else
                cMap(i,:) = [0.5, 0.5, 0.5];% gray
            end
        end

        % Define frequency bands and corresponding colors
        bands = {
            [0, 4],      [1 0 0];      % red for 0–4 Hz
            [5, 9],      [0 0 1];      % blue for 5–9 Hz
            [10, 15],    [0 1 0];      % green for 10–15 Hz
            [16, 20],    [0.5 0.5 0.5] % gray for 16–20 Hz
        };

        nBands = size(bands, 1);
        bandMeans = zeros(1, nBands);

        % Calculate the mean PSD (in dB) for each frequency band
        for b = 1:nBands
            range = bands{b, 1};
            idx = (freq >= range(1)) & (freq <= range(2));
            if any(idx)
                bandMeans(b) = mean(psd_dB(idx));
            else
                bandMeans(b) = NaN;  % or 0 if you prefer
            end
        end

        % 7) Create a bar chart for this bin's PSD
        %    "FaceColor=flat" means we can assign colors via CData.
        %    "FaceAlpha=0.3" for partial transparency.
        b = bar(G.A10, freq, psd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3, ...
                'DisplayName', sprintf('Epoch# %d', centerBin));

        % Assign the color map
        b.CData = cMap;
        % 9) Axis labeling
        xlabel(G.A10, 'Frequency (Hz)');
        xlim(G.A10, freqRange);
        grid(G.A10,'on');
        legend(G.A10,'show');  % show each bin in the legend

        hold(G.A10, 'off');



        % % -------
        % We'll track min/max across all bins to standardize the y-scale
        if centerBin <=5
            centerBin_MinusOne = 1;
            centerBin_MinusTwo = 1;
            centerBin_MinusThree = 1;
            centerBin_MinusFour = 1;
        end
        current_data_MO = getBinData(G, centerBin_MinusOne);


        [MOfreq, MOpsdVals] = computeBinPSD(current_data_MO, fs, windowSeconds, overlapPercent);

        MOmask = (MOfreq>= freqRange(1)) & (MOfreq <= freqRange(2));
        MOfreq = MOfreq(MOmask);
        MOpsdVals = MOpsdVals(MOmask);
        MOpsd_dB = MOpsdVals;
        printBandMeans(MOfreq, MOpsd_dB);
        MOpsd_dB = MOpsdVals;
        % 6) Color-code each frequency point in freq by 4 bands:
        %    0-4 Hz => red, 5-9 Hz => blue, 10-15 Hz => green, 16-20 Hz => gray
        cMap = zeros(length(MOfreq), 3);  % each row = RGB
        for i = 1:length(MOfreq)
            f = freq(i);
            if f >= 0  && f <= 4
                cMap(i,:) = [1, 0, 0];      % red
            elseif f > 4  && f <= 9
                cMap(i,:) = [0, 0, 1];      % blue
            elseif f > 9 && f <= 15
                cMap(i,:) = [0, 1, 0];      % green
            else
                cMap(i,:) = [0.5, 0.5, 0.5];% gray
            end
        end

        % Define frequency bands and corresponding colors
        bands = {
            [0, 4],      [1 0 0];      % red for 0–4 Hz
            [5, 9],      [0 0 1];      % blue for 5–9 Hz
            [10, 15],    [0 1 0];      % green for 10–15 Hz
            [16, 20],    [0.5 0.5 0.5] % gray for 16–20 Hz
        };

        nBands = size(bands, 1);
        bandMeans = zeros(1, nBands);

        % Calculate the mean PSD (in dB) for each frequency band
        for b = 1:nBands
            range = bands{b, 1};
            idx = (MOfreq >= range(1)) & (MOfreq <= range(2));
            if any(idx)
                bandMeans(b) = mean(MOpsd_dB(idx));
            else
                bandMeans(b) = NaN;  % or 0 if you prefer
            end
        end

        % 7) Create a bar chart for this bin's PSD
        %    "FaceColor=flat" means we can assign colors via CData.
        %    "FaceAlpha=0.3" for partial transparency.
        MOb = bar(G.A20, MOfreq, MOpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3, ...
                'DisplayName', sprintf('Epoch# %d', centerBin_MinusOne));

        % Assign the color map
        MOb.CData = cMap;
        % 9) Axis labeling
        xlim(G.A20, freqRange);
        grid(G.A20,'on');
        legend(G.A20,'show');  % show each bin in the legend
        ylim(G.A20, [0, meanVal + 3*stdDev]);
        hold(G.A20, 'off');

        % CurrentBin - 2
        current_data_MT = getBinData(G, centerBin_MinusTwo);

        [MTfreq, MTpsdVals] = computeBinPSD(current_data_MT, fs, windowSeconds, overlapPercent);

        MTmask = (MTfreq>= freqRange(1)) & (MTfreq <= freqRange(2));
        MTfreq = MTfreq(MTmask);
        MTpsdVals = MTpsdVals(MTmask);
        MTpsd_dB = MTpsdVals;
        printBandMeans(MTfreq, MTpsd_dB);
        MTpsd_dB = MTpsdVals;
        % 6) Color-code each frequency point in freq by 4 bands:
        %    0-4 Hz => red, 5-9 Hz => blue, 10-15 Hz => green, 16-20 Hz => gray
        cMap = zeros(length(MTfreq), 3);  % each row = RGB
        for i = 1:length(MTfreq)
            f = freq(i);
            if f >= 0  && f <= 4
                cMap(i,:) = [1, 0, 0];      % red
            elseif f > 4  && f <= 9
                cMap(i,:) = [0, 0, 1];      % blue
            elseif f > 9 && f <= 15
                cMap(i,:) = [0, 1, 0];      % green
            else
                cMap(i,:) = [0.5, 0.5, 0.5];% gray
            end
        end

        % Define frequency bands and corresponding colors
        bands = {
            [0, 4],      [1 0 0];      % red for 0–4 Hz
            [5, 9],      [0 0 1];      % blue for 5–9 Hz
            [10, 15],    [0 1 0];      % green for 10–15 Hz
            [16, 20],    [0.5 0.5 0.5] % gray for 16–20 Hz
        };

        nBands = size(bands, 1);
        bandMeans = zeros(1, nBands);

        % Calculate the mean PSD (in dB) for each frequency band
        for b = 1:nBands
            range = bands{b, 1};
            idx = (MTfreq >= range(1)) & (MTfreq <= range(2));
            if any(idx)
                bandMeans(b) = mean(MTpsd_dB(idx));
            else
                bandMeans(b) = NaN;  % or 0 if you prefer
            end
        end

        % 7) Create a bar chart for this bin's PSD
        %    "FaceColor=flat" means we can assign colors via CData.
        %    "FaceAlpha=0.3" for partial transparency.
        MTb = bar(G.A19, MTfreq, MTpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3, ...
                'DisplayName', sprintf('Epoch# %d', centerBin_MinusTwo));

        % Assign the color map
        MTb.CData = cMap;
        % 9) Axis labeling
        xlim(G.A19, freqRange);
        grid(G.A19,'on');
        legend(G.A19,'show');  % show each bin in the legend
        ylim(G.A19, [0, meanVal + 3*stdDev]);
        hold(G.A19, 'off');


          % CurrentBin - 3
        current_data_MTH = getBinData(G, centerBin_MinusThree);

        [MTHfreq, MTHpsdVals] = computeBinPSD(current_data_MTH, fs, windowSeconds, overlapPercent);

        MTHmask = (MTHfreq>= freqRange(1)) & (MTHfreq <= freqRange(2));
        MTHfreq = MTHfreq(MTHmask);
        MTHpsdVals = MTHpsdVals(MTHmask);
        MTHpsd_dB = MTHpsdVals;
        printBandMeans(MTHfreq, MTHpsd_dB);
        MTHpsd_dB = MTHpsdVals;
        % 6) Color-code each frequency point in freq by 4 bands:
        %    0-4 Hz => red, 5-9 Hz => blue, 10-15 Hz => green, 16-20 Hz => gray
        cMap = zeros(length(MTHfreq), 3);  % each row = RGB
        for i = 1:length(MTHfreq)
            f = freq(i);
            if f >= 0  && f <= 4
                cMap(i,:) = [1, 0, 0];      % red
            elseif f > 4  && f <= 9
                cMap(i,:) = [0, 0, 1];      % blue
            elseif f > 9 && f <= 15
                cMap(i,:) = [0, 1, 0];      % green
            else
                cMap(i,:) = [0.5, 0.5, 0.5];% gray
            end
        end

        % Define frequency bands and corresponding colors
        bands = {
            [0, 4],      [1 0 0];      % red for 0–4 Hz
            [5, 9],      [0 0 1];      % blue for 5–9 Hz
            [10, 15],    [0 1 0];      % green for 10–15 Hz
            [16, 20],    [0.5 0.5 0.5] % gray for 16–20 Hz
        };

        nBands = size(bands, 1);
        bandMeans = zeros(1, nBands);

        % Calculate the mean PSD (in dB) for each frequency band
        for b = 1:nBands
            range = bands{b, 1};
            idx = (MTHfreq >= range(1)) & (MTHfreq <= range(2));
            if any(idx)
                bandMeans(b) = mean(MTHpsd_dB(idx));
            else
                bandMeans(b) = NaN;  % or 0 if you prefer
            end
        end

        % 7) Create a bar chart for this bin's PSD
        %    "FaceColor=flat" means we can assign colors via CData.
        %    "FaceAlpha=0.3" for partial transparency.
        MTHb = bar(G.A18, MTHfreq, MTHpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3, ...
                'DisplayName', sprintf('Epoch# %d', centerBin_MinusThree));

        % Assign the color map
        MTHb.CData = cMap;
        % 9) Axis labeling
        xlim(G.A18, freqRange);
        grid(G.A18,'on');
        legend(G.A18,'show');  % show each bin in the legend
        ylim(G.A18, [0, meanVal + 3*stdDev]);
        hold(G.A18, 'off');


        % CurrentBin - 4
        current_data_MF = getBinData(G, centerBin_MinusFour);

        [MFfreq, MFpsdVals] = computeBinPSD(current_data_MF, fs, windowSeconds, overlapPercent);

        MFmask = (MFfreq>= freqRange(1)) & (MFfreq <= freqRange(2));
        MFfreq = MFfreq(MFmask);
        MFpsdVals = MFpsdVals(MFmask);
        MFpsd_dB = MFpsdVals;
        printBandMeans(MFfreq, MFpsd_dB);
        MFpsd_dB = MFpsdVals;
        % 6) Color-code each frequency point in freq by 4 bands:
        %    0-4 Hz => red, 5-9 Hz => blue, 10-15 Hz => green, 16-20 Hz => gray
        cMap = zeros(length(MFfreq), 3);  % each row = RGB
        for i = 1:length(MFfreq)
            f = freq(i);
            if f >= 0  && f <= 4
                cMap(i,:) = [1, 0, 0];      % red
            elseif f > 4  && f <= 9
                cMap(i,:) = [0, 0, 1];      % blue
            elseif f > 9 && f <= 15
                cMap(i,:) = [0, 1, 0];      % green
            else
                cMap(i,:) = [0.5, 0.5, 0.5];% gray
            end
        end

        % Define frequency bands and corresponding colors
        bands = {
            [0, 4],      [1 0 0];      % red for 0–4 Hz
            [5, 9],      [0 0 1];      % blue for 5–9 Hz
            [10, 15],    [0 1 0];      % green for 10–15 Hz
            [16, 20],    [0.5 0.5 0.5] % gray for 16–20 Hz
        };

        nBands = size(bands, 1);
        bandMeans = zeros(1, nBands);

        % Calculate the mean PSD (in dB) for each frequency band
        for b = 1:nBands
            range = bands{b, 1};
            idx = (MFfreq >= range(1)) & (MFfreq <= range(2));
            if any(idx)
                bandMeans(b) = mean(MFpsd_dB(idx));
            else
                bandMeans(b) = NaN;  % or 0 if you prefer
            end
        end

        % 7) Create a bar chart for this bin's PSD
        %    "FaceColor=flat" means we can assign colors via CData.
        %    "FaceAlpha=0.3" for partial transparency.
        MFb = bar(G.A17, MFfreq, MFpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3, ...
                'DisplayName', sprintf('Epoch# %d', centerBin_MinusFour));

        % Assign the color map
        MFb.CData = cMap;
        % 9) Axis labeling
        xlim(G.A17, freqRange);
        grid(G.A17,'on');
        legend(G.A17,'show');  % show each bin in the legend
        ylim(G.A17, [0, meanVal + 3*stdDev]);
        hold(G.A17, 'off');

        % CurrentBin + 1
        current_data_PO = getBinData(G, centerBin_PlusOne);

        [POfreq, POpsdVals] = computeBinPSD(current_data_PO, fs, windowSeconds, overlapPercent);

        POmask = (POfreq>= freqRange(1)) & (POfreq <= freqRange(2));
        POfreq = POfreq(POmask);
        POpsdVals = POpsdVals(POmask);
        POpsd_dB = POpsdVals;
        printBandMeans(POfreq, POpsd_dB);
        POpsd_dB = POpsdVals;
        % 6) Color-code each frequency point in freq by 4 bands:
        %    0-4 Hz => red, 5-9 Hz => blue, 10-15 Hz => green, 16-20 Hz => gray
        cMap = zeros(length(POfreq), 3);  % each row = RGB
        for i = 1:length(POfreq)
            f = freq(i);
            if f >= 0  && f <= 4
                cMap(i,:) = [1, 0, 0];      % red
            elseif f > 4  && f <= 9
                cMap(i,:) = [0, 0, 1];      % blue
            elseif f > 9 && f <= 15
                cMap(i,:) = [0, 1, 0];      % green
            else
                cMap(i,:) = [0.5, 0.5, 0.5];% gray
            end
        end

        % Define frequency bands and corresponding colors
        bands = {
            [0, 4],      [1 0 0];      % red for 0–4 Hz
            [5, 9],      [0 0 1];      % blue for 5–9 Hz
            [10, 15],    [0 1 0];      % green for 10–15 Hz
            [16, 20],    [0.5 0.5 0.5] % gray for 16–20 Hz
        };

        nBands = size(bands, 1);
        bandMeans = zeros(1, nBands);

        % Calculate the mean PSD (in dB) for each frequency band
        for b = 1:nBands
            range = bands{b, 1};
            idx = (POfreq >= range(1)) & (POfreq <= range(2));
            if any(idx)
                bandMeans(b) = mean(POpsd_dB(idx));
            else
                bandMeans(b) = NaN;  % or 0 if you prefer
            end
        end

        % 7) Create a bar chart for this bin's PSD
        %    "FaceColor=flat" means we can assign colors via CData.
        %    "FaceAlpha=0.3" for partial transparency.
        POb = bar(G.A21, POfreq, POpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3, ...
                'DisplayName', sprintf('Epoch# %d', centerBin_PlusOne));

        % Assign the color map
        POb.CData = cMap;
        % 9) Axis labeling
        xlim(G.A21, freqRange);
        grid(G.A21,'on');
        legend(G.A21,'show');  % show each bin in the legend
        ylim(G.A21, [0, meanVal + 3*stdDev]);
        hold(G.A21, 'off');

        % CurrentBin +2
        current_data_MF = getBinData(G, centerBin_MinusFour);

        [MFfreq, MFpsdVals] = computeBinPSD(current_data_MF, fs, windowSeconds, overlapPercent);

        MFmask = (MFfreq>= freqRange(1)) & (MFfreq <= freqRange(2));
        MFfreq = MFfreq(MFmask);
        MFpsdVals = MFpsdVals(MFmask);
        MFpsd_dB = MFpsdVals;
        printBandMeans(MFfreq, MFpsd_dB);
        MFpsd_dB = MFpsdVals;
        % 6) Color-code each frequency point in freq by 4 bands:
        %    0-4 Hz => red, 5-9 Hz => blue, 10-15 Hz => green, 16-20 Hz => gray
        cMap = zeros(length(MFfreq), 3);  % each row = RGB
        for i = 1:length(MFfreq)
            f = freq(i);
            if f >= 0  && f <= 4
                cMap(i,:) = [1, 0, 0];      % red
            elseif f > 4  && f <= 9
                cMap(i,:) = [0, 0, 1];      % blue
            elseif f > 9 && f <= 15
                cMap(i,:) = [0, 1, 0];      % green
            else
                cMap(i,:) = [0.5, 0.5, 0.5];% gray
            end
        end

        % Define frequency bands and corresponding colors
        bands = {
            [0, 4],      [1 0 0];      % red for 0–4 Hz
            [5, 9],      [0 0 1];      % blue for 5–9 Hz
            [10, 15],    [0 1 0];      % green for 10–15 Hz
            [16, 20],    [0.5 0.5 0.5] % gray for 16–20 Hz
        };

        nBands = size(bands, 1);
        bandMeans = zeros(1, nBands);

        % Calculate the mean PSD (in dB) for each frequency band
        for b = 1:nBands
            range = bands{b, 1};
            idx = (MFfreq >= range(1)) & (MFfreq <= range(2));
            if any(idx)
                bandMeans(b) = mean(MFpsd_dB(idx));
            else
                bandMeans(b) = NaN;  % or 0 if you prefer
            end
        end

        % 7) Create a bar chart for this bin's PSD
        %    "FaceColor=flat" means we can assign colors via CData.
        %    "FaceAlpha=0.3" for partial transparency.
        MFb = bar(G.A17, MFfreq, MFpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3, ...
                'DisplayName', sprintf('Epoch# %d', centerBin_MinusFour));

        % Assign the color map
        MFb.CData = cMap;
        % 9) Axis labeling
        xlim(G.A17, freqRange);
        grid(G.A17,'on');
        legend(G.A17,'show');  % show each bin in the legend
        ylim(G.A17, [0, meanVal + 3*stdDev]);
        hold(G.A17, 'off');

        % CurrentBin + 1
        current_data_PO = getBinData(G, centerBin_PlusOne);

        [POfreq, POpsdVals] = computeBinPSD(current_data_PO, fs, windowSeconds, overlapPercent);

        POmask = (POfreq>= freqRange(1)) & (POfreq <= freqRange(2));
        POfreq = POfreq(POmask);
        POpsdVals = POpsdVals(POmask);
        POpsd_dB = POpsdVals;
        printBandMeans(POfreq, POpsd_dB);
        POpsd_dB = POpsdVals;
        % 6) Color-code each frequency point in freq by 4 bands:
        %    0-4 Hz => red, 5-9 Hz => blue, 10-15 Hz => green, 16-20 Hz => gray
        cMap = zeros(length(POfreq), 3);  % each row = RGB
        for i = 1:length(POfreq)
            f = freq(i);
            if f >= 0  && f <= 4
                cMap(i,:) = [1, 0, 0];      % red
            elseif f > 4  && f <= 9
                cMap(i,:) = [0, 0, 1];      % blue
            elseif f > 9 && f <= 15
                cMap(i,:) = [0, 1, 0];      % green
            else
                cMap(i,:) = [0.5, 0.5, 0.5];% gray
            end
        end

        % Define frequency bands and corresponding colors
        bands = {
            [0, 4],      [1 0 0];      % red for 0–4 Hz
            [5, 9],      [0 0 1];      % blue for 5–9 Hz
            [10, 15],    [0 1 0];      % green for 10–15 Hz
            [16, 20],    [0.5 0.5 0.5] % gray for 16–20 Hz
        };

        nBands = size(bands, 1);
        bandMeans = zeros(1, nBands);

        % Calculate the mean PSD (in dB) for each frequency band
        for b = 1:nBands
            range = bands{b, 1};
            idx = (POfreq >= range(1)) & (POfreq <= range(2));
            if any(idx)
                bandMeans(b) = mean(POpsd_dB(idx));
            else
                bandMeans(b) = NaN;  % or 0 if you prefer
            end
        end

        % 7) Create a bar chart for this bin's PSD
        %    "FaceColor=flat" means we can assign colors via CData.
        %    "FaceAlpha=0.3" for partial transparency.
        POb = bar(G.A21, POfreq, POpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3, ...
                'DisplayName', sprintf('Epoch# %d', centerBin_PlusOne));

        % Assign the color map
        POb.CData = cMap;
        % 9) Axis labeling
        xlim(G.A21, freqRange);
        grid(G.A21,'on');
        legend(G.A21,'show');  % show each bin in the legend
        ylim(G.A21, [0, meanVal + 3*stdDev]);
        hold(G.A21, 'off');

                % CurrentBin +2
        current_data_PT = getBinData(G, centerBin_PlusTwo);

        [PTfreq, PTpsdVals] = computeBinPSD(current_data_PT, fs, windowSeconds, overlapPercent);

        PTmask = (PTfreq>= freqRange(1)) & (PTfreq <= freqRange(2));
        PTfreq = PTfreq(PTmask);
        PTpsdVals = PTpsdVals(PTmask);
        PTpsd_dB = PTpsdVals;
        printBandMeans(MFfreq, PTpsd_dB);
        PTpsd_dB = PTpsdVals;
        % 6) Color-code each frequency point in freq by 4 bands:
        %    0-4 Hz => red, 5-9 Hz => blue, 10-15 Hz => green, 16-20 Hz => gray
        cMap = zeros(length(PTfreq), 3);  % each row = RGB
        for i = 1:length(PTfreq)
            f = freq(i);
            if f >= 0  && f <= 4
                cMap(i,:) = [1, 0, 0];      % red
            elseif f > 4  && f <= 9
                cMap(i,:) = [0, 0, 1];      % blue
            elseif f > 9 && f <= 15
                cMap(i,:) = [0, 1, 0];      % green
            else
                cMap(i,:) = [0.5, 0.5, 0.5];% gray
            end
        end

        % Define frequency bands and corresponding colors
        bands = {
            [0, 4],      [1 0 0];      % red for 0–4 Hz
            [5, 9],      [0 0 1];      % blue for 5–9 Hz
            [10, 15],    [0 1 0];      % green for 10–15 Hz
            [16, 20],    [0.5 0.5 0.5] % gray for 16–20 Hz
        };

        nBands = size(bands, 1);
        bandMeans = zeros(1, nBands);

        % Calculate the mean PSD (in dB) for each frequency band
        for b = 1:nBands
            range = bands{b, 1};
            idx = (PTfreq >= range(1)) & (PTfreq <= range(2));
            if any(idx)
                bandMeans(b) = mean(PTpsd_dB(idx));
            else
                bandMeans(b) = NaN;  % or 0 if you prefer
            end
        end

        % 7) Create a bar chart for this bin's PSD
        %    "FaceColor=flat" means we can assign colors via CData.
        %    "FaceAlpha=0.3" for partial transparency.
        PTb = bar(G.A22, PTfreq, PTpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3, ...
                'DisplayName', sprintf('Epoch# %d', centerBin_PlusTwo));

        % Assign the color map
        PTb.CData = cMap;
        % 9) Axis labeling
        xlim(G.A22, freqRange);
        grid(G.A22,'on');
        legend(G.A22,'show');  % show each bin in the legend
        ylim(G.A22, [0, meanVal + 3*stdDev]);
        hold(G.A22, 'off');

        % CurrentBin + 3
        current_data_PTH = getBinData(G, centerBin_PlusThree);

        [PTHfreq, PTHpsdVals] = computeBinPSD(current_data_PTH, fs, windowSeconds, overlapPercent);

        PTHmask = (PTHfreq>= freqRange(1)) & (PTHfreq <= freqRange(2));
        PTHfreq = PTHfreq(PTHmask);
        PTHpsdVals = PTHpsdVals(PTHmask);
        PTHpsd_dB = PTHpsdVals;
        printBandMeans(PTHfreq, PTHpsd_dB);
        PTHpsd_dB = PTHpsdVals;
        % 6) Color-code each frequency point in freq by 4 bands:
        %    0-4 Hz => red, 5-9 Hz => blue, 10-15 Hz => green, 16-20 Hz => gray
        cMap = zeros(length(PTHfreq), 3);  % each row = RGB
        for i = 1:length(PTHfreq)
            f = freq(i);
            if f >= 0  && f <= 4
                cMap(i,:) = [1, 0, 0];      % red
            elseif f > 4  && f <= 9
                cMap(i,:) = [0, 0, 1];      % blue
            elseif f > 9 && f <= 15
                cMap(i,:) = [0, 1, 0];      % green
            else
                cMap(i,:) = [0.5, 0.5, 0.5];% gray
            end
        end

        % Define frequency bands and corresponding colors
        bands = {
            [0, 4],      [1 0 0];      % red for 0–4 Hz
            [5, 9],      [0 0 1];      % blue for 5–9 Hz
            [10, 15],    [0 1 0];      % green for 10–15 Hz
            [16, 20],    [0.5 0.5 0.5] % gray for 16–20 Hz
        };

        nBands = size(bands, 1);
        bandMeans = zeros(1, nBands);

        % Calculate the mean PSD (in dB) for each frequency band
        for b = 1:nBands
            range = bands{b, 1};
            idx = (POfreq >= range(1)) & (POfreq <= range(2));
            if any(idx)
                bandMeans(b) = mean(PTHpsd_dB(idx));
            else
                bandMeans(b) = NaN;  % or 0 if you prefer
            end
        end

        % 7) Create a bar chart for this bin's PSD
        %    "FaceColor=flat" means we can assign colors via CData.
        %    "FaceAlpha=0.3" for partial transparency.
        PTHb = bar(G.A23, PTHfreq, PTHpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3, ...
                'DisplayName', sprintf('Epoch# %d', centerBin_PlusThree));

        % Assign the color map
        PTHb.CData = cMap;
        % 9) Axis labeling
        xlim(G.A23, freqRange);
        grid(G.A23,'on');
        legend(G.A23,'show');  % show each bin in the legend
        ylim(G.A23, [0, meanVal + 3*stdDev]);
        hold(G.A23, 'off');


                % CurrentBin + 4
        current_data_PF = getBinData(G, centerBin_PlusFour);

        [PFfreq, PFpsdVals] = computeBinPSD(current_data_PF, fs, windowSeconds, overlapPercent);

        PFmask = (PFfreq>= freqRange(1)) & (PFfreq <= freqRange(2));
        PFfreq = PFfreq(PFmask);
        PFpsdVals = PFpsdVals(PFmask);
        PFpsd_dB = PFpsdVals;
        printBandMeans(PFfreq, PFpsd_dB);
        PFpsd_dB = PFpsdVals;
        % 6) Color-code each frequency point in freq by 4 bands:
        %    0-4 Hz => red, 5-9 Hz => blue, 10-15 Hz => green, 16-20 Hz => gray
        cMap = zeros(length(PFfreq), 3);  % each row = RGB
        for i = 1:length(PFfreq)
            f = freq(i);
            if f >= 0  && f <= 4
                cMap(i,:) = [1, 0, 0];      % red
            elseif f > 4  && f <= 9
                cMap(i,:) = [0, 0, 1];      % blue
            elseif f > 9 && f <= 15
                cMap(i,:) = [0, 1, 0];      % green
            else
                cMap(i,:) = [0.5, 0.5, 0.5];% gray
            end
        end

        % Define frequency bands and corresponding colors
        bands = {
            [0, 4],      [1 0 0];      % red for 0–4 Hz
            [5, 9],      [0 0 1];      % blue for 5–9 Hz
            [10, 15],    [0 1 0];      % green for 10–15 Hz
            [16, 20],    [0.5 0.5 0.5] % gray for 16–20 Hz
        };

        nBands = size(bands, 1);
        bandMeans = zeros(1, nBands);

        % Calculate the mean PSD (in dB) for each frequency band
        for b = 1:nBands
            range = bands{b, 1};
            idx = (PFfreq >= range(1)) & (PFfreq <= range(2));
            if any(idx)
                bandMeans(b) = mean(PFpsd_dB(idx));
            else
                bandMeans(b) = NaN;  % or 0 if you prefer
            end
        end

        % 7) Create a bar chart for this bin's PSD
        %    "FaceColor=flat" means we can assign colors via CData.
        %    "FaceAlpha=0.3" for partial transparency.
        PFb = bar(G.A24, PFfreq, PFpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3, ...
                'DisplayName', sprintf('Epoch# %d', centerBin_PlusFour));

        % Assign the color map
        PFb.CData = cMap;
        % 9) Axis labeling
        xlim(G.A24, freqRange);
        grid(G.A24,'on');
        legend(G.A24,'show');  % show each bin in the legend
        ylim(G.A24, [0, meanVal + 3*stdDev]);
        hold(G.A24, 'off');

        % for binIdx = binStart : binEnd
        %     % 2) Extract data for this bin
        %     data = getBinData(G, binIdx);
        %
        %
        %
        %     if isempty(data)
        %         continue;  % skip invalid or out-of-bounds bin
        %     end
        %
        %     % 3) Compute PSD (freq, psdVals) via Welch
        %     [freq, psdVals] = computeBinPSD(data, fs, windowSeconds, overlapPercent);
        %
        %     localMax = max(psdVals);
        %     localMin = min(psdVals);
        %
        %
        %     % 4) Filter freq range
        %     mask = (freq >= freqRange(1)) & (freq <= freqRange(2));
        %     freq    = freq(mask);
        %     psdVals = psdVals(mask);
        %
        %     % 5) Convert to dB if desired
        %     % psd_dB = 10 * log10(psdVals);
        %     psd_dB = psdVals;
        %     printBandMeans(freq, psd_dB);
        %     psd_dB = psdVals;
        %
        %     % 6) Color-code each frequency point in freq by 4 bands:
        %     %    0-4 Hz => red, 5-9 Hz => blue, 10-15 Hz => green, 16-20 Hz => gray
        %     cMap = zeros(length(freq), 3);  % each row = RGB
        %     for i = 1:length(freq)
        %         f = freq(i);
        %         if f >= 0  && f <= 4
        %             cMap(i,:) = [1, 0, 0];      % red
        %         elseif f > 4  && f <= 9
        %             cMap(i,:) = [0, 0, 1];      % blue
        %         elseif f > 9 && f <= 15
        %             cMap(i,:) = [0, 1, 0];      % green
        %         else
        %             cMap(i,:) = [0.5, 0.5, 0.5];% gray
        %         end
        %     end
        %
        %     % Define frequency bands and corresponding colors
        %     bands = {
        %         [0, 4],      [1 0 0];      % red for 0–4 Hz
        %         [5, 9],      [0 0 1];      % blue for 5–9 Hz
        %         [10, 15],    [0 1 0];      % green for 10–15 Hz
        %         [16, 20],    [0.5 0.5 0.5] % gray for 16–20 Hz
        %     };
        %
        %     nBands = size(bands, 1);
        %     bandMeans = zeros(1, nBands);
        %
        %     % Calculate the mean PSD (in dB) for each frequency band
        %     for b = 1:nBands
        %         range = bands{b, 1};
        %         idx = (freq >= range(1)) & (freq <= range(2));
        %         if any(idx)
        %             bandMeans(b) = mean(psd_dB(idx));
        %         else
        %             bandMeans(b) = NaN;  % or 0 if you prefer
        %         end
        %     end
        %
        %     % 7) Create a bar chart for this bin's PSD
        %     %    "FaceColor=flat" means we can assign colors via CData.
        %     %    "FaceAlpha=0.3" for partial transparency.
        %     b = bar(G.A10, freq, psd_dB, 1, ...
        %             'EdgeColor','none', ...
        %             'FaceColor','flat', ...
        %             'FaceAlpha',0.3, ...
        %             'DisplayName', sprintf('Epoch# %d', binIdx));
        %
        %     % Assign the color map
        %     b.CData = cMap;
        %
        % end
        %
        % % 9) Axis labeling
        % xlabel(G.A10, 'Frequency (Hz)');
        % xlim(G.A10, freqRange);
        % grid(G.A10,'on');
        % legend(G.A10,'show');  % show each bin in the legend
        %
        % % 10) Standardize y-lim across all bins
        %     yRange = localMax - localMin;
        %     ymin   = localMin - 0.1*yRange;
        %     % ymax   = max(totalpsdVals);
        %     ymax = localMax;
        %     ylim(G.A10, [ymin, ymax]);
        %
        %
        % hold(G.A10, 'off');
    end

    function data = getBinData(G, binIdx)
        data = [];
        fs = 1 / G.dt;  % sampling freq
        epochSamples = round(G.epochLen * fs);

        % If binIdx < 1 or binIdx > G.nbins, might be out of range
        if binIdx < 1 || binIdx > G.nbins
            return; % just return empty
        end

        % Start sample
        startSample = (binIdx - 1)*epochSamples + 1;
        % End sample
        endSample   = binIdx*epochSamples;

        % clamp
        if startSample < 1
            startSample = 1;
        end
        if endSample > length(G.EEG)
            endSample = length(G.EEG);
        end

        data = G.EEG(startSample:endSample);
    end
    function [freq, psdVals] = computeBinPSD(data, fs, windowSec, overlapPercent)
        N = length(data);

        % sub-window length in samples:
        wlen = round(windowSec * fs);
        if wlen > N
            wlen = N;
        end

        noverlap = round(wlen * overlapPercent);
        if noverlap >= wlen
            noverlap = wlen-1;
        end

        [psdVals, freq] = pwelch(data, hann(wlen), noverlap, wlen, fs);
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
        if ~isfield(G, 'keyPressCount')
            G.keyPressCount = 0;
        end
        % Increment the counter
        G.keyPressCount = G.keyPressCount + 1;
        % Check if it's time to save
        if G.keyPressCount == 30
            saveFile();
            G.keyPressCount = 0;
        end
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
                % axes(G.A1);
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

            case {'r*','4','numpad4'} % set to REM*
                G.labels(G.index) = 4;
                G.unsavedChanges = 1;
                updateState;
                updatePlots;

                advance;

            case {'w*','5','numpad5'} % set to wake*
                G.labels(G.index) = 5;
                G.unsavedChanges = 1;
                updateState;
                updatePlots;

                advance;
            case {'s*','6','numpad6'} % set to NREM*
                G.labels(G.index) = 6;
                G.unsavedChanges = 1;
                updateState;
                updatePlots;

                advance;




            case 'f' % save file
                saveFile();

            case 'h' % show help menu
                % showHelp(G.A1, []);

            case 'multiply' % apply label to range of timepoints
                t = text(G.A5,18.75,-.9,sprintf(['Move the ROI\n',...
                    'boundaries,\nand then\ndouble-click it',...
                    '\nor press\nescape']),'Color','r');
                % axes(G.A1);
                % set(G.A1,'Clipping','off')
                % xl = xlim(G.A1);
                d = diff(xl);
                % roi = imrect(G.A1,[xl(1)+d/2 - d/24,-4.0287,d/12,4.9080]);
                rectPosition = round(wait(roi)./(G.epochLen/3600));
                roi.delete();
                % set(G.A1,'Clipping','on')
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
    function recenterEEG()
        % recenterEEG Re-center the EEG y-axis (G.A6a) so that the limits are symmetric about 0.
        %
        % This function checks the current y-axis limits (in volts) on G.A6a.
        % It then calculates the smaller of the absolute values of the lower and
        % upper limits and sets new limits as [-clipVal, clipVal].

        % Retrieve current y-axis limits from G.A6a (in volts)
        currentYLim = get(G.A6, 'YLim');

        % Compute the absolute values of the lower and upper limits
        lowerAbs = abs(currentYLim(1));
        upperAbs = abs(currentYLim(2));

        % Choose the smaller absolute value as the clipping value
        clipVal = min(lowerAbs, upperAbs);

        % Define new symmetric y-axis limits
        newYLim = [-clipVal, clipVal];

        % Apply the new y-axis limits to G.A6a and force manual mode
        set(G.A6, 'YLim', newYLim, 'YLimMode', 'manual');
        set(G.A6a, 'YLim', newYLim, 'YLimMode', 'manual');

    end

    function recenterEMG()
        % recenterEEG Re-center the EEG y-axis (G.A6a) so that the limits are symmetric about 0.
        %
        % This function checks the current y-axis limits (in volts) on G.A6a.
        % It then calculates the smaller of the absolute values of the lower and
        % upper limits and sets new limits as [-clipVal, clipVal].

        % Retrieve current y-axis limits from G.A6a (in volts)
        currentYLim = get(G.A7, 'YLim');

        % Compute the absolute values of the lower and upper limits
        lowerAbs = abs(currentYLim(1));
        upperAbs = abs(currentYLim(2));

        % Choose the smaller absolute value as the clipping value
        clipVal = min(lowerAbs, upperAbs);

        % Define new symmetric y-axis limits
        newYLim = [-clipVal, clipVal];

        % Apply the new y-axis limits to G.A6a and force manual mode
        set(G.A7, 'YLim', newYLim, 'YLimMode', 'manual');
        set(G.A7a, 'YLim', newYLim, 'YLimMode', 'manual');

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
        % G.emgYlim = [G.emgYlim(1)+.05*diff(G.emgYlim), G.emgYlim(2)-.05*diff(G.emgYlim)];
        % updatePlots;
        % defocus(src);
        currentYLim = get(G.A6, 'YLim');
        rangeY = diff(currentYLim);
        % Reduce the half-range by 15%
        newYLim = [currentYLim(1) + 0.15 * rangeY, currentYLim(2) - 0.15 * rangeY];
        set(G.A6, 'YLim', newYLim, 'YLimMode', 'manual');
        set(G.A6a, 'YLim', newYLim, 'YLimMode', 'manual');
        set(G.A6b, 'YLim', newYLim, 'YLimMode', 'manual');
        % Re-center the axis symmetrically about 0
        recenterEEG;
        updatePlots;
        defocus(src);
end

function fct_zoomoutEEG(src, ~)
        % G.emgYlim = [G.emgYlim(1)-.05*diff(G.emgYlim), G.emgYlim(2)+.05*diff(G.emgYlim)];
        % updatePlots;
        % defocus(src);
        currentYLim = get(G.A6, 'YLim');
        rangeY = diff(currentYLim);
        % Reduce the half-range by 15%
        newYLim = [currentYLim(1) - 0.15 * rangeY, currentYLim(2) + 0.15 * rangeY];
        set(G.A6, 'YLim', newYLim, 'YLimMode', 'manual');
        set(G.A6a, 'YLim', newYLim, 'YLimMode', 'manual');
        set(G.A6b, 'YLim', newYLim, 'YLimMode', 'manual');

        % Re-center the axis symmetrically about 0
        recenterEEG;
        updatePlots;
        defocus(src);
end

    function fct_zoominEMG(src,~)

        % G.emgYlim = [G.emgYlim(1)+.05*diff(G.emgYlim), G.emgYlim(2)-.05*diff(G.emgYlim)];
        % updatePlots;
        % defocus(src);
        currentYLim = get(G.A7, 'YLim');
        rangeY = diff(currentYLim);
        % Reduce the half-range by 15%
        newYLim = [currentYLim(1) + 0.15 * rangeY, currentYLim(2) - 0.15 * rangeY];
        set(G.A7, 'YLim', newYLim, 'YLimMode', 'manual');
        set(G.A7a, 'YLim', newYLim, 'YLimMode', 'manual');
        set(G.A7b, 'YLim', newYLim, 'YLimMode', 'manual');
        % Re-center the axis symmetrically about 0
        recenterEMG;
        updatePlots;
        defocus(src);
    end

    function fct_zoomoutEMG(src,~)
        %
        % G.emgYlim = [G.emgYlim(1)-.05*diff(G.emgYlim), G.emgYlim(2)+.05*diff(G.emgYlim)];
        % updatePlots;
        % defocus(src);
        currentYLim = get(G.A7, 'YLim');
        rangeY = diff(currentYLim);
        % Reduce the half-range by 15%
        newYLim = [currentYLim(1) - 0.15 * rangeY, currentYLim(2) + 0.15 * rangeY];
        set(G.A7, 'YLim', newYLim, 'YLimMode', 'manual');
        set(G.A7a, 'YLim', newYLim, 'YLimMode', 'manual');
        set(G.A7b, 'YLim', newYLim, 'YLimMode', 'manual');
        % Re-center the axis symmetrically about 0
        recenterEEG;
        updatePlots;
        defocus(src);
    end

    %% Helper functions

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
        if ~isempty(setdiff(unique(x.labels), 1:6)) % in the range 1:4
            msgbox('Error: labels variable must be in the range 1:6')
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
