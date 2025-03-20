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

% VARIABLES/CONSTANTS
TARGET_SAMPLE_RATE = 256;

% EEG LOWPASS VALUE
EEG_LOW_PASS_CUTOFF = 25;

% EMG BANDPASS VALUES
EMG_LOW_BANDPASS_CUTOFF = 25;
EMG_HIGH_BANDPASS_CUTOFF = 50;

%EEG & EMG Y-LIM Multiplier
YLIM_MULTIPLIER = 2.2;

% SPECTROGRAM - SHOW VALUES UP TO 21Hz
SPECTROGRAM_CUTOFF = 21;

% Welch analysis parameters
windowSeconds  = 4;      % window length for Welch
overlapPercent = 0.9;    % 90% overlap
freqRange      = [0 20]; % show 0-20 Hz

% Default # Bins to show (15-1)/2 = 7 bins on each side [+ current bin]
DEFAULT_SHOW = 15;



% PSD - DELTA - RED, THETA - BLUE, ALPHA - GREEN, ELSE - GRAY
DELTA_COLOR = [1, 0, 0];
THETA_COLOR = [0, 0, 1];
ALPHA_COLOR = [0, 1, 0];
ELSE_COLOR = [0.5, 0.5, 0.5];

PADDING = 0.02;

LABEL_STAGE_ALL_COLOR = [1 1 1; .47 .67 .19; .14 .2 .57; 0.996 0.758 0.039; 0.25 0.45 0.05; 0.05 0.15 0.50;0.85 0.65 0.03];
WHITE = LABEL_STAGE_ALL_COLOR(1, :); %BACKGROUND
GREEN = LABEL_STAGE_ALL_COLOR(2, :); %REM
BLUE = LABEL_STAGE_ALL_COLOR(3, :); %WAKE
YELLOW = LABEL_STAGE_ALL_COLOR(4, :); %NREM
DARK_GREEN = LABEL_STAGE_ALL_COLOR(5, :); %REM-Artifact
DARK_BLUE = LABEL_STAGE_ALL_COLOR(6, :); %WAKE_Artifact
DARK_YELLOW = LABEL_STAGE_ALL_COLOR(7, :); %NREM_Artifact

% Compute Target Sample Rate, and epoch Length
G.originalSR = SR; % EEG/EMG sampling rate
G.SR = TARGET_SAMPLE_RATE; % sampling rate used when calculating spectrogram and processed EMG
G.epochLen = epochLen; % length of one epoch (spectrogram column) in seconds

if length(EEG) ~= length(EMG)
    message = 'ERROR: EEG and EMG are different lengths';
    return
end

% Filter EEG, EMG data - Low pass -> EEG, bandpass -> EMG
EEG = lowpass(EEG, EEG_LOW_PASS_CUTOFF, SR);
EMG = bandpass(EMG,[EMG_LOW_BANDPASS_CUTOFF EMG_HIGH_BANDPASS_CUTOFF], SR);

G.EEG = EEG - mean(EEG);
G.EMG = EMG - mean(EMG);
clear('EMG','EEG');


% create spectrogram and process EMG using G.SR == TARGET_SAMPLE_RATE
[spec, tAxis, fAxis] = createSpectrogram(standardizeSR(G.EEG, G.originalSR, G.SR), G.SR, G.epochLen);
G.processedEMG = processEMG(standardizeSR(G.EMG, G.originalSR, G.SR), G.SR, G.epochLen);
G.cappedEMG = G.processedEMG;
emgCap = mean(G.cappedEMG) + 3*std(G.cappedEMG);
G.cappedEMG(G.cappedEMG > emgCap) = emgCap;

% set various parameters
G.show = DEFAULT_SHOW; %  default number of bins to display on screen
G.dt = 1/G.originalSR; % duration of each EEG/EMG sample in seconds
G.advance = 1; % whether to advance automatically when a state is assigned
G.colors = [WHITE;GREEN;BLUE;YELLOW;DARK_GREEN;DARK_BLUE;DARK_YELLOW]; %colors for sleep stages
G.mid = ceil(G.show/2); % important for plotting the current time marker - middle of G.show
G.savepath = ''; % where to save the sleep stage labels
G.nbins = length(tAxis); % total number of time bins in the recording
G.unsavedChanges = 0; % whether anything has changed since user last saved
fs = 1 / G.dt;  % sampling rate

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
showFreqs = find(fAxis <= SPECTROGRAM_CUTOFF); % only show frequencies under 21 Hz
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
G.eegYlim = [-YLIM_MULTIPLIER*yl, YLIM_MULTIPLIER*yl];
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
    'Position', [0, 0, 1, 1],'KeyPressFcn',@keypress,...
    'Menubar', 'none','Color', 'w', 'Name', 'AccuSleep_viewer');


% Spectrogram
G.A3 = axes('Units', 'Normalized', 'Position', [0.05 0.95 0.87 0.045]);
G.A3.Toolbar.Visible = 'off';
set(gca, 'FontSize', 10, 'LineWidth', 2, 'XTick', [], 'YTick', []);
ylabel(G.A3, 'Spec.');

% Upper time point indicator
G.A8 = axes('Units', 'Normalized', 'Position', [0.05 0.935  0.87 0.015],'XTick',[],'YTick',[]);
G.A8.Toolbar.Visible = 'off';

% Upper sleep stage labels
G.A1 = axes('Units', 'Normalized', 'Position', [0.05 0.88 0.87 0.05]);
G.A1.Toolbar.Visible = 'off';

% EEG signal - 1 Minute
G.A6a = axes('Units', 'Normalized', 'Position', [0.05 0.775, 0.87 .09]);
G.A6a.Toolbar.Visible = 'off';
ylabel(G.A6a, 'EEG-1');
set(G.A6a,'XTick',[])

% % EMG signal - 1 Minute
G.A7a = axes('Units', 'Normalized', 'Position', [0.05 0.67 0.87 .09]);
G.A7a.Toolbar.Visible = 'off';
ylabel(G.A7a, 'EMG-1');
set(G.A7a,'XTick',[])

% EEG signal  - main
G.A6 = axes('Units', 'Normalized', 'Position', [0.05 0.565 0.87 .09]);
G.A6.Toolbar.Visible = 'off';
ylabel(G.A6, 'EEG');

% EMG signal  - main
G.A7 = axes('Units', 'Normalized', 'Position', [0.05 0.445 0.87 .09]);
G.A7.Toolbar.Visible = 'off';
ylabel(G.A7, 'EMG');
set(G.A7,'XTick',[])

% EEG signal + 1 minute
G.A6b = axes('Units', 'Normalized', 'Position', [0.05 0.34 0.87 .09]);
G.A6b.Toolbar.Visible = 'off';
ylabel(G.A6b, 'EEG+1');
set(G.A6b,'XTick',[])

% EMG signal + 1 minute
G.A7b = axes('Units', 'Normalized', 'Position', [0.05 0.235 0.87 .09]);
G.A7b.Toolbar.Visible = 'off';
ylabel(G.A7b, 'EMG+1');
set(G.A7b,'XTick',[])

% -4 PSD
G.A17 = axes('Units', 'Normalized', 'Position', [0.05 0.13 0.085 0.09]);
G.A17.Toolbar.Visible = 'off';
% -3 PSD
G.A18 = axes('Units', 'Normalized', 'Position', [0.145 0.13 0.085 0.09]);
G.A18.Toolbar.Visible = 'off';
% -2 PSD
G.A19 = axes('Units', 'Normalized', 'Position', [0.24 0.13 0.085 0.09]);
G.A19.Toolbar.Visible = 'off';
% -1 PSD
G.A20 = axes('Units', 'Normalized', 'Position', [0.34 0.13 0.085 0.09]);
G.A20.Toolbar.Visible = 'off';
% Current PSD - 0 PSD
G.A10 = axes('Units', 'Normalized', 'Position', [0.44 0.13 0.085 0.09]);
G.A10.Toolbar.Visible = 'off';
% +1 PSD
G.A21 = axes('Units', 'Normalized', 'Position', [0.54 0.13 0.085 0.09]);
G.A21.Toolbar.Visible = 'off';
% +2 PSD
G.A22 = axes('Units', 'Normalized', 'Position', [0.645 0.13 0.085 0.09]);
G.A22.Toolbar.Visible = 'off';
% +3 PSD
G.A23 = axes('Units', 'Normalized', 'Position', [0.74 0.13 0.085 0.09]);
G.A23.Toolbar.Visible = 'off';
% +4 PSD
G.A24 = axes('Units', 'Normalized', 'Position', [0.84 0.13 0.085 0.09]);
G.A24.Toolbar.Visible = 'off';

% Lower sleep stage labels
G.A9 = axes('Units', 'Normalized', 'Position', [0.05 0.01 0.87 0.09]);


% Plot the EEG data
globalMin_EEG = min(G.EEG);
globalMax_EEG = max(G.EEG);
globalRange_EEG = globalMax_EEG - globalMin_EEG;
padding_EEG = PADDING * globalRange_EEG;
stableMin_EEG = globalMin_EEG - padding_EEG;
stableMax_EEG = globalMax_EEG + padding_EEG;
% --- In updatePlots, apply the stable EEG axis limits ---
set(G.A6, 'YLim', [stableMin_EEG, stableMax_EEG], 'YLimMode', 'manual');
set(G.A6a, 'YLim', [stableMin_EEG, stableMax_EEG], 'YLimMode', 'manual');
set(G.A6b, 'YLim', [stableMin_EEG, stableMax_EEG], 'YLimMode', 'manual');

% Plot the EMG data
globalMin_EMG = min(G.EMG);
globalMax_EMG = max(G.EMG);
globalRange_EMG = globalMax_EMG - globalMin_EMG;
padding_EMG = PADDING * globalRange_EMG;
stableMin_EMG = globalMin_EMG - padding_EMG;
stableMax_EMG = globalMax_EMG + padding_EMG;
% --- In updatePlots, apply the stable EMG axis limits ---
set(G.A7, 'YLim', [stableMin_EMG, stableMax_EMG], 'YLimMode', 'manual');
set(G.A7a, 'YLim', [stableMin_EMG, stableMax_EMG], 'YLimMode', 'manual');
set(G.A7b, 'YLim', [stableMin_EMG, stableMax_EMG], 'YLimMode', 'manual');

linkaxes([G.A1,G.A3, G.A8], 'x'); % upper panel x axes should stay linked

% Buttons
G.helpbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized','BackgroundColor',[1 .8 .8],...
    'Position',[.93 .01 .062 .04],'String','Help','Callback',@showHelp,'FontSize',9,...
    'ToolTip','Show help menu (H)');
G.savebtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized','BackgroundColor',[.8 1 .8],...
    'Position',[.93 .055 .062 .04],'String','Save labels','Callback',@saveCallback,'FontSize',9,...
    'ToolTip','Save labels to file (F)');
G.loadbtn = uicontrol(WIN,'Style','pushbutton', 'Units','normalized',...
    'Position',[.93 .10 .062 .03],'String','Load labels','Callback',@loadFile,'FontSize',9,...
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
% Text Boxes
G.sumDelta = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [.925 .18 .07 .035], 'String', '', ...
    'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
G.sumTheta = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [.925 .16 .07 .035], 'String', '', ...
    'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
G.thetaRatio = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [.925 .14 .07 .035], 'String', '', ...
    'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
G.emgProcessed = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [.925 .32 .07 .035], 'String', '', ...
    'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
% keep track of the current timepoint
% G.eegIntegral = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
%     'Position', [.925 .32 .07 .035], 'String', '', ...
%     'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
% G.emgIntegral = uicontrol(WIN, 'Style', 'text', 'Units', 'normalized', ...
%     'Position', [.925 .32 .07 .035], 'String', '', ...
%     'FontSize', 9, 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
% keep track of the current timepoint

G.index = 1; % index of current time point
G.timepointS = G.specTs(G.index); % current time point, in seconds
G.timepointH = G.specTh(G.index); % current time point, in hours

% plot spectrogram
clim(G.A3,G.caxis1);
imagesc(G.A3,G.specTh, fAxis(showFreqs), G.spectrogram',G.caxis1);
axis(G.A3, 'xy')
colormap(G.A3,G.colormap);
G.lims = xlim(G.A3);
set(G.A3, 'YTick', 0:5:20, 'YTickLabel', 0:5:20);
ylabel(G.A3, 'Spec.');
set(G.A3, 'XTick', [], 'XTickLabel', [], 'XColor', 'none');

updateState;

% Compute EMG Integral
% EMG_integral = computeSignalIntegral(G.processedEMG,G.SR);
% G.EMG_Integral = EMG_integral;


[~, totalpsdVals] = computeBinPSD(G.EEG, G.SR, 4, 0.9);
PSD_Scale = max(totalpsdVals) * 1.35;
ylim(G.A20, [0, PSD_Scale]);

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

        set(G.A9, 'XTickLabel', [],'XTick',[], 'YTick', [1 2 3], 'YTickLabel', {'REM', 'WAKE', 'NREM'});

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
        % Get the EMG Integral of the current bin
        % curEMGIntegral = G.EMG_Integral(ii);
        % curEMGProcessed = G.cappedEMG(ii);
        y_min_EMG = min(curEMG);
        y_max_EMG = max(curEMG);
        yr_EMG = y_max_EMG - y_min_EMG;
        padding_EMG = PADDING * yr_EMG;

        % Set the EMG integral to the text box
        % set(G.emgIntegral, 'String', sprintf('EMG Integral: %.2f ', curEMGIntegral));

        line(G.A7, t, curEMG, 'Color', 'k', 'LineWidth', 1);

        currentYLim = get(G.A7, 'YLim');  % currentYLim = [bottom, top]
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
        padding_EEG = PADDING * yr_EEG;
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

            % Plot the previous EEG data (in volts)
            line(G.A6a, t_prev, curEEG_prev, 'Color','k', 'LineWidth', 1);

            % EMG-previous minutes -1 minute
            cla(G.A7a);
            hold(G.A7a, 'on');
            xlim(G.A7a, [t_prev(1)-G.dt, t_prev(end)]);

            % Compute dynamic y-axis limits for previous EEG (in volts)
            curEMG_prev = G.EMG(ii_prev);
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

            % Plot the future EEG data (in volts)
            line(G.A6b, t_future, curEEG_future, 'Color','k', 'LineWidth', 1);


            cla(G.A7b);
            hold(G.A7b, 'on');
            xlim(G.A7b, [t_future(1)-G.dt, t_future(end)]);

            % Compute dynamic y-axis limits for the future EEG window (in volts)
            curEMG_future = G.EMG(ii_future);

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
        set(G.A1, 'XTickLabel', [],'XTick',[], 'YTick', [1 2 3], 'YTickLabel', {'REM', 'WAKE', 'NREM'});
    end
    % function integralSignal = computeSignalIntegral(signal, fs, options)
    %     % Default options
    %     if nargin < 3 || isempty(options)
    %         options.baselineCorrect = false;
    %         options.normalize = false;
    %         options.method = 'trapz';
    %     end
    %
    %     % Ensure signal is a column vector
    %     if isrow(signal)
    %         signal = signal';
    %     end
    %
    %     % Get number of samples
    %     nSamples = length(signal);
    %
    %     % Create time vector (in seconds)
    %     t = (0:nSamples-1) / fs;
    %
    %     % Baseline correction (remove mean if requested)
    %     if options.baselineCorrect
    %         signal = signal - mean(signal);
    %         fprintf('Signal baseline corrected (mean removed).\n');
    %     end
    %
    %     % Compute integral based on method
    %     switch lower(options.method)
    %         case 'trapz'
    %             % Use trapezoidal numerical integration (MATLAB's trapz)
    %             integralSignal = trapz(t, signal);
    %             % Reshape to match signal length for cumulative effect
    %             integralSignal = cumsum(signal) * (t(2) - t(1)); % Scale by time step
    %         case 'cumsum'
    %             % Use cumulative sum (simpler but less accurate for non-uniform data)
    %             dt = 1 / fs;  % Time step
    %             integralSignal = cumsum(signal) * dt;
    %         otherwise
    %             error('Unsupported integration method. Use ''trapz'' or ''cumsum''.');
    %     end
    %
    %     % Normalize by total time if requested
    %     if options.normalize
    %         totalTime = t(end) - t(1);
    %         integralSignal = integralSignal / totalTime;
    %         fprintf('Integral normalized by total time (%f seconds).\n', totalTime);
    %     end
    %
    %     % Ensure output is a column vector
    %     if isrow(integralSignal)
    %         integralSignal = integralSignal';
    %     end
    %
    %     fprintf('Integral computed with %d samples, fs = %d Hz.\n', nSamples, fs);
    % end
    % function integralSignal = computeSignalIntegral(signal, fs)
    %     % Ensure the signal is a column vector
    %     if isrow(signal)
    %         signal = signal';
    %     end
    %
    %     % Calculate the time step based on sampling frequency
    %     dt = 1 / fs;
    %
    %     % Compute the cumulative integral using the cumulative sum
    %     integralSignal = cumsum(signal) * dt;
    % end
    function CalculateBandSums(freq, psd_dB)
        % Define frequency bands and a corresponding label (for printing)
        bands = {
            [0, 4]
            [5, 9]
            [10, 15]
            [16, 20]
        };

        nBands = size(bands, 1);
        bandMeans = zeros(1, nBands);

        % Calculate the mean PSD (in dB) for each frequency band
        for b = 1:nBands
            range = bands{b, 1};
            idx = (freq >= range(1)) & (freq <= range(2));
            if any(idx)
                bandMeans(b) = sum(psd_dB(idx));
            else
                bandMeans(b) = NaN;
            end
        end


        % Calculate the mean power in the 0-4 Hz band
        sumPower_0_4Hz = sum(psd_dB(freq >= 0 & freq <= 4));
        sumPower_5_9Hz = sum(psd_dB(freq >=5 & freq <=9));
        tRatio = sumPower_5_9Hz/sumPower_0_4Hz;

        % Scientific notiation with .2e
        % Update the text box
        set(G.sumDelta, 'String', sprintf('Delta: %.2e ', sumPower_0_4Hz));
        set(G.sumTheta, 'String', sprintf('Theta: %.2e ', sumPower_5_9Hz));
        set(G.thetaRatio, 'String', sprintf('Theta Ratio: %.2f ', tRatio));

        % Print the mean power in the 0-4 Hz band
        fprintf('sum power in 0-4 Hz: %.10f dB\n', sumPower_0_4Hz);
        fprintf('sum power in 5-9 Hz: %.10f dB\n', sumPower_5_9Hz);
        fprintf('Theta-Ratio: %.10f dB\n', tRatio);
    end
    function showAllEpochsPSD(G)
        % 1) Clear G.A10 and hold on
        cla(G.A10);
        hold(G.A10, 'on');

        % PSD Bins
        centerBin = G.index;
        offsets = -4:4;            % Offsets from -4 to +4
        bins = centerBin + offsets;
        bins = max(bins, 1);         % Ensure no index is below 1

        % Now assign individual variables if needed
        centerBin_MinusFour   = bins(1);
        centerBin_MinusThree  = bins(2);
        centerBin_MinusTwo    = bins(3);
        centerBin_MinusOne    = bins(4);
        centerBin             = bins(5);
        centerBin_PlusOne     = bins(6);
        centerBin_PlusTwo     = bins(7);
        centerBin_PlusThree   = bins(8);
        centerBin_PlusFour    = bins(9);

        % Current Bin Data
        current_data = getBinData(G, centerBin);
        [freq, psdVals] = computeBinPSD(current_data, fs, windowSeconds, overlapPercent);
        mask = (freq>= freqRange(1)) & (freq <= freqRange(2));
        freq = freq(mask);
        psdVals = psdVals(mask);
        psd_dB = psdVals;
        CalculateBandSums(freq, psd_dB);

        % Color-code each frequency point in freq by 3 bands:
        %    (0<=4 Hz => red), (4<=9 Hz => blue), (9<=15 Hz => green),
        %    Else Gray
        cMap = zeros(length(freq), 3);  % each row = RGB
        for i = 1:length(freq)
            f = freq(i);
            if f >= 0  && f <= 4
                cMap(i,:) = DELTA_COLOR;      % red
            elseif f > 4  && f <= 9
                cMap(i,:) = THETA_COLOR;      % blue
            elseif f > 9 && f <= 15
                cMap(i,:) = ALPHA_COLOR;      % green
            else
                cMap(i,:) = ELSE_COLOR;       % gray
            end
        end

        b = bar(G.A10, freq, psd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3);
        b.CData = cMap;
        xlim(G.A10, freqRange);
        grid(G.A10,'on');
        ylim(G.A10, [0, PSD_Scale]);
        set(G.A10, 'YTickLabel', []);
        set(G.A10, 'XTickLabel', []);
        hold(G.A10, 'off');


        % CurrentBin - 1
        current_data_MO = getBinData(G, centerBin_MinusOne);
        [MOfreq, MOpsdVals] = computeBinPSD(current_data_MO, fs, windowSeconds, overlapPercent);
        MOmask = (MOfreq>= freqRange(1)) & (MOfreq <= freqRange(2));
        MOfreq = MOfreq(MOmask);
        MOpsdVals = MOpsdVals(MOmask);
        MOpsd_dB = MOpsdVals;
        MOb = bar(G.A20, MOfreq, MOpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3);
        MOb.CData = cMap;
        xlim(G.A20, freqRange);
        grid(G.A20,'on');
        ylim(G.A20, [0, PSD_Scale]);
        set(G.A20, 'YTickLabel', []);
        set(G.A20, 'XTickLabel', []);
        hold(G.A20, 'off');

        % CurrentBin - 2
        current_data_MT = getBinData(G, centerBin_MinusTwo);
        [MTfreq, MTpsdVals] = computeBinPSD(current_data_MT, fs, windowSeconds, overlapPercent);
        MTmask = (MTfreq>= freqRange(1)) & (MTfreq <= freqRange(2));
        MTfreq = MTfreq(MTmask);
        MTpsdVals = MTpsdVals(MTmask);
        MTpsd_dB = MTpsdVals;
        MTb = bar(G.A19, MTfreq, MTpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3);
        MTb.CData = cMap;
        xlim(G.A19, freqRange);
        grid(G.A19,'on');
        ylim(G.A19, [0, PSD_Scale]);
        set(G.A19, 'YTickLabel', []);
        set(G.A19, 'XTickLabel', []);
        hold(G.A19, 'off');

          % CurrentBin - 3
        current_data_MTH = getBinData(G, centerBin_MinusThree);
        [MTHfreq, MTHpsdVals] = computeBinPSD(current_data_MTH, fs, windowSeconds, overlapPercent);
        MTHmask = (MTHfreq>= freqRange(1)) & (MTHfreq <= freqRange(2));
        MTHfreq = MTHfreq(MTHmask);
        MTHpsdVals = MTHpsdVals(MTHmask);
        MTHpsd_dB = MTHpsdVals;
        MTHb = bar(G.A18, MTHfreq, MTHpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3);
        MTHb.CData = cMap;
        xlim(G.A18, freqRange);
        grid(G.A18,'on');
        ylim(G.A18, [0, PSD_Scale]);
        set(G.A18, 'YTickLabel', []);
        set(G.A18, 'XTickLabel', []);
        hold(G.A18, 'off');


        % CurrentBin - 4
        current_data_MF = getBinData(G, centerBin_MinusFour);
        [MFfreq, MFpsdVals] = computeBinPSD(current_data_MF, fs, windowSeconds, overlapPercent);
        MFmask = (MFfreq>= freqRange(1)) & (MFfreq <= freqRange(2));
        MFfreq = MFfreq(MFmask);
        MFpsdVals = MFpsdVals(MFmask);
        MFpsd_dB = MFpsdVals;
        MFb = bar(G.A17, MFfreq, MFpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3);
        MFb.CData = cMap;
        xlim(G.A17, freqRange);
        grid(G.A17,'on');
        ylim(G.A17, [0, PSD_Scale]);
        hold(G.A17, 'off');
        xlabel(G.A17, 'Frequency (Hz)');
        ylabel(G.A17, 'Power');


        % CurrentBin + 1
        current_data_PO = getBinData(G, centerBin_PlusOne);
        [POfreq, POpsdVals] = computeBinPSD(current_data_PO, fs, windowSeconds, overlapPercent);
        POmask = (POfreq>= freqRange(1)) & (POfreq <= freqRange(2));
        POfreq = POfreq(POmask);
        POpsdVals = POpsdVals(POmask);
        POpsd_dB = POpsdVals;
        POb = bar(G.A21, POfreq, POpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3);
        POb.CData = cMap;
        xlim(G.A21, freqRange);
        grid(G.A21,'on');
        ylim(G.A21, [0, PSD_Scale]);
        set(G.A21, 'YTickLabel', []);
        set(G.A21, 'XTickLabel', []);
        hold(G.A21, 'off');



        % CurrentBin + 2
        current_data_PT = getBinData(G, centerBin_PlusTwo);
        [PTfreq, PTpsdVals] = computeBinPSD(current_data_PT, fs, windowSeconds, overlapPercent);
        PTmask = (PTfreq>= freqRange(1)) & (PTfreq <= freqRange(2));
        PTfreq = PTfreq(PTmask);
        PTpsdVals = PTpsdVals(PTmask);
        PTpsd_dB = PTpsdVals;
        PTb = bar(G.A22, PTfreq, PTpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3);
        PTb.CData = cMap;
        xlim(G.A22, freqRange);
        grid(G.A22,'on');
        ylim(G.A22, [0, PSD_Scale]);
        set(G.A22, 'YTickLabel', []);
        set(G.A22, 'XTickLabel', []);
        hold(G.A22, 'off');

        % CurrentBin + 3
        current_data_PTH = getBinData(G, centerBin_PlusThree);
        [PTHfreq, PTHpsdVals] = computeBinPSD(current_data_PTH, fs, windowSeconds, overlapPercent);
        PTHmask = (PTHfreq>= freqRange(1)) & (PTHfreq <= freqRange(2));
        PTHfreq = PTHfreq(PTHmask);
        PTHpsdVals = PTHpsdVals(PTHmask);
        PTHpsd_dB = PTHpsdVals;
        PTHb = bar(G.A23, PTHfreq, PTHpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3);
        PTHb.CData = cMap;
        xlim(G.A23, freqRange);
        grid(G.A23,'on');
        ylim(G.A23, [0, PSD_Scale]);
        set(G.A23, 'YTickLabel', []);
        set(G.A23, 'XTickLabel', []);
        hold(G.A23, 'off');


        % CurrentBin + 4
        current_data_PF = getBinData(G, centerBin_PlusFour);
        [PFfreq, PFpsdVals] = computeBinPSD(current_data_PF, fs, windowSeconds, overlapPercent);
        PFmask = (PFfreq>= freqRange(1)) & (PFfreq <= freqRange(2));
        PFfreq = PFfreq(PFmask);
        PFpsdVals = PFpsdVals(PFmask);
        PFpsd_dB = PFpsdVals;
        PFb = bar(G.A24, PFfreq, PFpsd_dB, 1, ...
                'EdgeColor','none', ...
                'FaceColor','flat', ...
                'FaceAlpha',0.3);
        PFb.CData = cMap;
        xlim(G.A24, freqRange);
        grid(G.A24,'on');
        ylim(G.A24, [0, PSD_Scale]);
        set(G.A24, 'YTickLabel', []);
        set(G.A24, 'XTickLabel', []);
        hold(G.A24, 'off');
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
        % Check if it's time to save file - 30 minutes is 450 keystrokes
        % 4 seconds x 15 = 60 [1 minute]
        % 15 * 30 = 450 keystrokes
        if G.keyPressCount == 450
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
                showHelp(G.A1, []);

            case 'multiply' % apply label to range of timepoints
                t = text(G.A5,18.75,-.9,sprintf(['Move the ROI\n',...
                    'boundaries,\nand then\ndouble-click it',...
                    '\nor press\nescape']),'Color','r');
                rectPosition = round(wait(roi)./(G.epochLen/3600));
                roi.delete();
                if isempty(rectPosition)
                    delete(t);
                    return
                end

                [label,~] = listdlg('PromptString','Set label to:',...
                    'SelectionMode','single',...
                    'ListString',{'REM', 'WAKE','NREM','Undefined'});
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
        clim(G.A3,G.caxis1)
        updatePlots;
        defocus(src);
    end

    function dimSpect(src,~)

        G.cmax = G.cmax + G.cmax/10;
        G.caxis1 = [G.caxis1(1), G.cmax];
        clim(G.A3,G.caxis1)
        updatePlots;
        defocus(src);
    end

    function recenterEEG()
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

    function fct_showmenu(src, ~)

       options = [1, 3, 5, 7, 9, 15];
        G.show = options(src.Value);
        G.mid = ceil(G.show/2);
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
