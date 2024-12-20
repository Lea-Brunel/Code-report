% Sleep scoring analysis
clc
clear all

% ---add path with Matlab functions---
p = genpath('D:\matlab');
addpath(p)

% ---add path with the data file---
path = 'C:\Users\tofb6\OneDrive\Documents\Stage INCI 2024\m34\060924 closed economy 24 hours small ROI';
cd(path)

load('data.mat')

SleepScoringGUI
%% plot hypnogram

load('hypnogram.mat');

figure
tiledlayout(3,1)

% First plot
ax1 = nexttile;
plot(data.signal(4,:));
xlim([0 size(data.signal, 2)])
ylim([-0.05 0.05])

% Second plot
ax2 = nexttile;
plot(data.signal(5,:));
xlim([0 samples])
ylim([-0.05 0.05])

% Third plot
ax3 = nexttile;
imagesc(hypnogram);

linkaxes([ax1 ax2 ax3],'x')

%% Plot hypnogram and feeding bouts

figure
imagesc(hypnogram);
hold on
plot(data.fed, 'r')
hold off

%% Get sleep scoring 1h before and after initiation of meal 

state = [];
state.hypnogram = hypnogram;
Fs_bin = 2.5;
hour = 60*60*Fs_bin; 

for i = 1:numel(batch_event)
    % Indices before the event
    start_before = max(1, batch_event(i) - hour); % Prevent indices before start of table
    end_before = batch_event(i);

    % Indices after the event
    start_after = batch_event(i);
    end_after = min(length(state.hypnogram), batch_event(i) + hour); % Prevent indices beyond table end

    % Initialize rows with NaN
    state.before(i, :) = NaN(1, hour + 1);
    state.after(i, :) = NaN(1, hour + 1);

    % Fill before with actual data
    actual_before = state.hypnogram(start_before:end_before);
    valid_start = hour + 1 - length(actual_before) + 1; % Adjust for shorter valid data
    state.before(i, valid_start:end) = actual_before;

    % Fill after with actual data
    actual_after = state.hypnogram(start_after:end_after);
    state.after(i, 1:length(actual_after)) = actual_after;
end

%% Get proportion of wake before and after meal

for i = 1:numel(batch_event)
    state.before_wake_prop(i) = nnz(state.before(i,:)==1)/sum(~isnan(state.before(i,:))); % nnz=sum
    state.after_wake_prop(i) = nnz(state.after(i,:)==1)/ sum(~isnan(state.after(i,:))); 
end
    
%% Separate event with sleep/without sleep 10 min before meal initiation 

time = 10 ; % in minutes

% get sleep scoring 10 min before meal
for i = 1:numel(batch_event)
    state.tenmin(i,:) = state.before(i,end-time*60*Fs_bin:end);
end

% separate meals
a=1;
b=1;

for i = 1:numel(batch_event)
    if all(state.tenmin(i,:)==1) % 1 = awake
        for channel=1:size(data.good_sig_batch,1)
            state.nosleep10(channel,a,:) = data.good_sig_batch(channel,i,:);
        end
        state.batchnosleep10(a,:) = data.fed_batch(i);
        a = a+1;
    else % 2 and 3 = sleep
        for channel=1:size(data.good_sig_batch,1)
            state.withsleep10(channel,b,:) = data.good_sig_batch(channel,i,:);
        end
        state.batchsleep10(b,:) = data.fed_batch(i);
        b = b+1;
    end
end

%% plot mean signal sleep/no sleep

F = Fs_bin*60; % frames for 1 minute (150)
range = 10;
x = (1:2*range*F+1)/F;
centre = (-range*F : range*F)/F;

ch = 8;

figure
subplot(211)

% without sleep
input = squeeze(state.nosleep10(ch,:,:));
input1 = mean(input(:,:),1,'omitnan');
input2 = std(input(:,:),[],1,'omitnan')/sqrt(size(input,1));
shadedErrorBar(centre,input1,input2,'lineProps','-b','transparent',1);
xline(0,'--r')
ylabel('$\Delta F/F_0$', 'Interpreter', 'latex')
title('Mean signal no sleep before')

subplot(212)

% with sleep
input = squeeze(state.withsleep10(ch,:,:));
input1 = mean(input(:,:),1,'omitnan');
input2 = std(input(:,:),[],1,'omitnan')/sqrt(size(input,1));
shadedErrorBar(centre,input1,input2,'lineProps','-b','transparent',1);

xline(0,'--r')
title('Mean signal with sleep before')
xlabel('minutes before feeding')
ylabel('$\Delta F/F_0$', 'Interpreter', 'latex')

%% Test relationship

for channel = 1:size(state.nosleep10,1)
    input = squeeze(state.nosleep10(channel,:,:));
    input1 = (input(:,1500) - input(:,1200))/300 ; % rate of change
    data.rise(channel, :, :) = input1;
    mdl = fitglm(state.nosleep10(channel,:), closed.fed_batch, 'linear', 'Distribution', 'poisson');
    
    disp(mdl)
end
