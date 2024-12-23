%% Data processing and analysis
% ------- READ DATA -------
clc
clear all

% ---add path with Matlab functions---
p = genpath('D:\matlab');
addpath(p)

% ---add path with the data file---
path = 'C:\Users\tofb6\OneDrive\Documents\Stage INCI 2024\m34\110924 free feeding 24 hours';
cd(path)
%% Process data
fdata = lvm_import('test_FED_m1');
temp = fdata.Segment1.data(:,2);

data = lvm_import('test_GREEN_m1');

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

if round(numel(data)/2)<numel(m_even)
    temp = [temp; nan(numel(m_even)*2-numel(temp),1)];
end
fed = temp(m_even);

if numel(m_even)< numel(m_odd)
    signal = [signal nan(25,1)];
    fed = [fed; nan];
elseif numel(m_even)> numel(m_odd)
    background = [background nan(25,1)];
end

% ---save extracted raw data ---
%save('data_raw','signal','background')

%%
% ------- BIN AND NORMALIZE -------

bin_step = 4; % bin means regrouping consecutive points to reduce noise, here one bin contains 4 pts
Fs= 10; % frequency of acquisition in Hz (number of frame per sec)
Fs_bin = Fs/bin_step; % new frequency after binning
percentile_norm=10; % percentile to normalize to

data = []; % data will contain every values of intermediates calculations and final signal to group everything in the same place

% open bin_norm_data function
% bin_norm_data( data , channels , bin step , type of normalization (1=10th percentile kind with deltaF/F), percentile to normalize to)
data.s_norm = bin_norm_data(signal', (1:24), bin_step, 1, percentile_norm)';
data.b_norm = bin_norm_data(background', (1:24), bin_step, 1, percentile_norm)'; % normalizing sig and backgrd put them on the same scale to allow substraction (it's the equivalent of the multiplication by the coef in the other method)
data.fed = bin_norm_data(fed, (1), bin_step, 0)'; % just bin the data without normalizing (no need) to have the same frequency and nbr of frames

% ------- SUBSTACT BACKGROUND -------

data.signal = data.s_norm - data.b_norm;

%% ISOLATE BEHAVIOR SIGNALS
% ------- ISOLATE EVENTS -------

k = 1;
event=[]; % will contain the time point of every feeding bout
for count = 2:numel(data.fed)

    if and(data.fed(count)>1, data.fed(count-1)<1) % clean the signal bc sometimes if the nose poke is long it will count as 2 consecutive pts when it's actually one
        % if the present pts is >1 (nose poke) and the previous one was <1 (no nose poke) then count
        event(k)= count;
        k = k+1;
    end
end

% ------- IDENTIFY MEALS -------
% group of feeding-bouts that happen in a short time window : 1 minute

batch_event(1)=event(1); %store the first value
k = 2;
for count = 2:numel(event)
    if event(count-1)<event(count)-Fs_bin*60
        batch_event(k)= event(count);
        k = k+1;
    end
end

% ------- CUT SIGNAL AROUND MEALS (batch) -------

data.sig_batch=[];
range = 10;
for count = 1:numel(batch_event)

    if batch_event(count)-Fs_bin*60*range<0 %if batch - 10 minutes <0 : checks if the window around the event extends before the start of the signal
        diff_range = abs(batch_event(count)-Fs_bin*60*range)+1; % count the number of missing values (it's negative from the start so abs takes the absolute value) and +1 bc we are indexing at -1
        data.sig_batch(:,count,:) = [nan(24,diff_range) data.signal(:,1:batch_event(count)+Fs_bin*60*range)]; % adds NaN were theres no signal
    else
        data.sig_batch(:,count,:) = data.signal(:,batch_event(count)-Fs_bin*60*range:batch_event(count)+Fs_bin*60*range);
    end
end

%% Get meal information : size and intervalle

data.fed_batch = []; % nbr of feeding-bout in a batch
data.batch_range = []; % tims btw two batch
data.batch_range(1) = nan;

for count = 2:numel(batch_event)

    % ------- TIME PASSED SINCE LAST FEEDING -------
    data.batch_range(count) = (batch_event(count) - batch_event(count-1));

    % ------- CALCULATE THE NBR OF COLLECTED PELLETS -------
    data.fed_batch(count-1) = sum(and( event > batch_event(count-1), event < batch_event(count))) +1;
end

% to include number of consumed pellets in the last batch_event
data.fed_batch(count) = sum( and(event> batch_event(count), event< numel(data.fed)))+1;

%% Save all data

%save('data','data')
load("data.mat")

%% Store channel correspondance (brain region)
num_channels = size(data.s_norm,1);
channel_names = cell(1, num_channels);

% loop for the channel's name
for i = 1:num_channels
    channel_names{i} = ['Channel ' num2str(i)];
end

% store name of the brain regions
channel_correspondance = {'Left LO Lateral Orbital Cortex', 'Left M2 Secondary Motor Cortex', ...
    'Left VO Ventral Orbital Cortex', 'Left M2 Secondary Motor Cortex',     'Left PrL Prelimbic Cortex',...
    'Left Cg1 Cingulate Cortex Area 1','Right Cg1 Cingulate Cortex Area 1','Right PrL Prelimbic Cortex',...
    'Right M2 Secondary Motor Cortex',    'Right VO Ventral Orbital Cortex','Right M2 Secondary Motor Cortex', ...
    'Right LO Lateral Orbital Cortex','Left PO Posterior Thalamic Nuclear Group', ...
    'Left GrDG Granule Cell Layer of the Dentate Gyrus','Left CA1 Field CA1 of the Hippocampus', ...
    'Left CL Centrolateral Thalamic Nucleus',    'Left MD Mediodorsal Thalamic Nucleus', ...
    'Left RSGc Retrosplenial Granular Cortex c region','Right RSGc Retrosplenial Granular Cortex c region', ...
    'Right MD Mediodorsal Thalamic Nucleus','Right CL Centrolateral Thalamic Nucleus','Right CA1 Field CA1 of the Hippocampus', ...
    'Right GrDG Granule Cell Layer of the Dentate Gyrus','Right PO Posterior Thalamic Nuclear Group'};

% store short name of brain regions
channel_correspondance_abr = {'Left LO', 'Left M2', 'Left VO', 'Left M2', 'Left PrL', ...
    'Left Cg1','Right Cg1','Right PrL','Right M2', 'Right VO','Right M2','Right LO', ...
    'Left PO','Left GrDG','Left CA1', 'Left CL', 'Left MD','Left RSGc','Right RSGc', ...
    'Right MD','Right CL','Right CA1', 'Right GrDG','Right PO'};

% create dictionnary
channel_dict = containers.Map(channel_names, channel_correspondance); % to call a name : channel_dict('Channel n')
channel_dict_abr = containers.Map(channel_names, channel_correspondance_abr);

%% Plot
% ------- SIGNAL VS BACKGROUND ONE CHANNEL -------

Fp= 60*60*Fs_bin; % 1 hour
t = 1/Fp:1/Fp:numel(data.s_norm(1,:))/Fp; % x axis (time)
ch=8;

figure
plot(t,data.s_norm(ch,:)) % signal

hold on
plot(t,data.b_norm(ch,:),'-g') % background

xlim([min(t) max(t)])
xlabel('Time (h)')
ylabel('$\Delta F/F_0$', 'Interpreter','latex')
legend('470nm', '405nm')
% legend('470 nm signal', '405 nm signal (isosbestic)')
hold off

ax = gca; % Obtenir l'axe courant  
ax.XColor = 'k'; % keep x line visible (black)
ax.YColor = 'k'; % keep y line visible  
ax.Box = 'off';  % remove the box around the plot
title(channel_correspondance_abr(ch))

%%
% ------- SIGNAL AND FEEDING BOUT 24h FOR ONE CHANNEL -------

Fp= 60*60*Fs_bin; % 1 hour
t = 1/Fp:1/Fp:numel(data.s_norm(1,:))/Fp;
ch = 8;

% creating binary array to plot feeding events 
for i = 1:size(data.fed,2)
    if data.fed(1,i)>1
        val = 1;
    else
        val = 0;
    end
    fedclean(i) = val;
end

figure
tiledlayout(2,1)

% Plot signal
ax1 = nexttile;
plot(t,data.signal(ch,:), Color='b')
xlim([min(t) max(t)])
%xlabel('Time (h)')
ylabel('$\Delta F/F_0$', 'Interpreter', 'latex'); 
title(ax1, channel_correspondance_abr(ch))
ax = gca; 
ax.XColor = 'k'; 
ax.YColor = 'k'; 
ax.Box = 'off';  

% Plot feeding bouts
ax2 = nexttile;
plot(t,fedclean, Color='k')
xlim([min(t) max(t)])
ylim(ax2, [-0.05 1])
xlabel('Time (h)')
%title(ax2, channel_correspondance_abr(8))
ax = gca; % Obtenir l'axe courant  
ax.XColor = 'k'; % keep x line visible (black)
ax.YColor = 'k'; % keep y line visible  
ax.Box = 'off';  % remove the box around the plot
ax.TickDir = 'out'; % put ticks belox the ax
yticks([0,1])
yticklabels({'Not feeding', 'Feeding'})

linkaxes([ax1 ax2],'x')
%%
% ------- AVERAGE SIGNAL AROUND MEAL FOR ALL CHANNELS -------

% open shadedErrorBar function

F = Fs_bin*60; % frames for 1 minute (150)
range = 10; % time interval
x = (1:2*range*F+1)/F; % x axis
centre = (-range*F : range*F)/F; % to centre the x axis on the feeding onset (0 = feeding)

rows=6 ; cols=4;
figure
for j = 1:(rows*cols)
    subplot(rows,cols,j)
    input = squeeze(data.sig_batch(j,:,:));
    input1 = mean(input(:,:),1,'omitnan');
    input2 = std(input(:,:),[],1,'omitnan')/sqrt(size(input,1));
    shadedErrorBar(centre,input1,input2,'lineProps','-b','transparent',1); 
    xlim([min(centre) max(centre)]);
    %xlabel('Time (minutes)')
    %ylabel('ΔF/F0')
    xlabel('Time from feeding onset (min)');
    ylabel('$\Delta F/F_0$', 'Interpreter', 'latex');  % latex is the name of scientific notation, must be between $
    hold on
    xline(0 ,'--r'); % dashed line at feeding pts
    hold off
    title(channel_correspondance_abr(j)) % channel name / to add channel nbr ['Channel ' num2str(j)], 
end
%sgtitle(sprintf('fluorescence mean across events | Intervals before and after (minutes) : %d', range))

%% Select channels with signal
% Create dictionnary with only the good channels (with a signal)
selected_channels = [1:12, 18:21, 23, 24];
good_channel_correspondance = cell(1,size(selected_channels,2)); % to store the name of the channels

for channel = 1:size(selected_channels,2)
    good_channel_correspondance{channel} = ['Channel ' num2str(selected_channels(channel))] ;
end 

channel_names_cut = cell(1,size(selected_channels,2));
for channel = 1:size(selected_channels,2)
    channel_names_cut{channel} = channel_correspondance_abr{selected_channels(channel)} ;
end 

%% Plot
% ------- PLOT MEAN SIGNAL AROUND EVENTS FOR GOOD CHANNELS ONLY -------

% store signal for the good channels
selected_channels = [1:12, 18:21, 23, 24];
data.good_sig_batch = data.sig_batch(selected_channels,:,:);

F = Fs_bin*60; % frames for 1 minute
range = 10;
centre = (-range*F : range*F)/F; 

rows=6 ; cols=3;
figure
for j = 1:(rows*cols)
    subplot(rows,cols,j)
    input = squeeze(data.good_sig_batch(j,:,:));
    input1 = mean(input(:,:),1,'omitnan');
    input2 = std(input(:,:),[],1,'omitnan')/sqrt(size(input,1));
    shadedErrorBar(centre,input1,input2,'lineProps','-b','transparent',1); 
    xlim([min(centre) max(centre)]);
    %xlabel('Time (minutes)')
    %ylabel('ΔF/F0')
    xlabel('Time from feeding onset (min)');
    ylabel('$\Delta F/F_0$', 'Interpreter', 'latex');   
    hold on
    % 
    % input = squeeze(closed.good_sig_batch(j,:,:));
    % input1 = mean(input(:,:),1,'omitnan');
    % input2 = std(input(:,:),[],1,'omitnan')/sqrt(size(input,1));
    % shadedErrorBar(centre,input1,input2,'lineProps','-g','transparent',1); 

    xline(0 ,'--r'); 
    hold off
    title(channel_names_cut(j)) % channel name / to add channel nbr ['Channel ' num2str(j)], 
end
% legend('free feeding', 'closed economy')

%%
% ------- HEATMAP GOOD CHANNELS -------
mu = squeeze(mean(data.good_sig_batch,2));
figure
imagesc(mu);
colorbar;
xline(size(mu,2)/2, '--r');
yticks(1:18);
yticklabels(channel_names_cut);
title('Heatmap of the signals around start of the feeding state for each channel')

%% Store rate of change for good channels

data.rise = [];

for channel = 1:size(data.good_sig_batch,1)
    input = squeeze(data.good_sig_batch(channel,:,:));
    input1 = (input(:,1500) - input(:,1150))/350 ;
    data.rise(channel, :, :) = input1;
end

%% Poisson regression for nbr of pellets x rate

for channel = 1:size(data.good_sig_batch,1)
    mdl = fitglm(data.rise(channel,:), data.fed_batch, 'linear', 'Distribution', 'poisson'); % modele poisson regression
    D = mdl.Deviance; % residual deviance

    mdl_null = fitglm(data.rise(channel,:), data.fed_batch, 'constant', 'Distribution', 'poisson'); % modele null
    D0 = mdl_null.Deviance; % null deviance

    r = 1-mdl.Deviance/mdl_null.Deviance; % pseudo R² qui mesure la qualité de l'ajustement du modele
    disp(mdl)
end

%% plot poisson regression
ch = 8; % chose channel to plot
mdl = fitglm(data.rise(ch,:), data.fed_batch, 'linear', 'Distribution', 'poisson'); 

x_pred = linspace(min(data.rise(ch,:)), max(data.rise(ch,:)), 100)';
% y_pred = predict(mdl, x_pred); % prediction pts accoring to poisson coefficients

[predicted_counts, CI] = predict(mdl, x_pred); % to get confidence interval bounds and predicted pts

% Plot observed vs predicted pts + confidence bounds
figure
scatter(data.rise(ch,:), data.fed_batch, 'x'); % Observed points

hold on;
% plot(x_pred, y_pred, 'r', 'LineWidth', 2, 'DisplayName', 'Poisson regression'); % Regression line
plot(x_pred, predicted_counts, 'r'); % regression line
plot(x_pred, CI(:,1), 'k--'); % Lower bound
plot(x_pred, CI(:,2), 'k--'); % Upper bound

xlabel('Number of pellets retrieved')
ylabel('Rate of change in fluorescence')
title('Poisson regression fit')
legend('Observed Data', 'Regression Line', '95% conf. bounds')
hold off;

%% linear regression time since last meal x rise

for channel = 1:size(data.good_sig_batch, 1)
    mdl = fitlm(data.batch_range, data.rise(channel,:));
    disp(mdl)
end

%% plot for final figure of correlation

F = Fs_bin*60; 
range = 10;
centre = (-range*F : range*F)/F; 
ch=8;

figure
% ------- Mean signal -------
subplot(3,1,1)

input = squeeze(data.sig_batch(ch,:,:));
input1 = mean(input(:,:),1,'omitnan');
input2 = std(input(:,:),[],1,'omitnan')/sqrt(size(input,1));

shadedErrorBar(centre,input1,input2,'lineProps','-b','transparent',1);
xlim([min(centre) max(centre)]);
xlabel('Time from feeding onset (minutes)')
ylabel('$\Delta F/F_0$', 'Interpreter', 'latex'); 
hold on
xline(0, '--r');
hold off
title(['Mean fluorescence signal aroud meal in the ' channel_names_cut(ch)])

% ------- Poisson regression -------
subplot(3,1,2)

mdl = fitglm(data.rise(ch,:), data.fed_batch, 'linear', 'Distribution', 'poisson'); 

x_pred = linspace(min(data.rise(ch,:)), max(data.rise(ch,:)), 100)'; 
[predicted_counts, CI] = predict(mdl, x_pred);

scatter(data.rise(ch,:), data.fed_batch, 'x'); % Observed points
hold on;
plot(x_pred, predicted_counts, 'r'); % Regression pts
plot(x_pred, CI(:,1), 'k-.'); % Lower bound
plot(x_pred, CI(:,2), 'k-.'); % Upper bound

ax = gca;
ax.Box = 'on';

ylabel('Number of pellets retrieved')
xlabel('Rate of change in fluorescence')
title('Poisson regression fit')
legend('Observed Data', 'Regression Line', '95% conf. bounds')
hold off;

% ------- Linear model -------
subplot(3,1,3)

range_batch = data.batch_range./F; % to have in minutes and not in frame nbr
mdl2 = fitlm(range_batch,data.rise(ch,:));
plt = plot(mdl2);
plt(3).LineStyle = "-."; % confidence bound line
ylabel('Rate of change in fluorescence')
xlabel('Time since the last meal (minutes)')

%disp(mdl2)
title('Linear fit to time passed since the last feeding')

%% Histogram
% ------- NUMBER OF PELLETS PER MEAL -------
figure
histogram(data.fed_batch, 'EdgeColor','w', 'FaceColor',"#f9a762", 'NumBins',9, 'BinLimits',[1 10])
xticks(1:10)
ylim([0 6.5])
%yticks([0, 10, 20, 30, 40, 50]) % for closed ec
xlabel('Number of pellets per meal')
% ylabel('Occurence')

% ------- TIME SINCE LAST MEAL -------
figure
histogram((data.batch_range./F), 'EdgeColor', 'w', 'FaceColor','#c09292','NumBins',11, 'BinLimits',[0 120])
yticks(1:4)
ylim([0 4.5])
xlabel('Time since last meal (minutes)')

%% mean + SEM of pellets and interval

input = closed.fed_batch; % or batch_range/F;

moy = mean(input, 'omitnan')
SEM = std(input, 'omitnan')/sqrt(numel(input))
