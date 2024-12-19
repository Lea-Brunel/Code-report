function y = bin_norm_data(x, limit_channels , bin_size, if_norm, p)

% args:
% x is input fluorescence 1 or 2D array samples x channels
% limit_channels limit calculation to given number of channels
% bin_size in samples
% normalize data: df/f to percentile if 1 , soft MIN-MAX if 2, 
% moving average non-normalized if 3
% p percentile to normalize to 

[k,l] = size(x);
count_channel=1;

% --- choose channels from PFC --- 
% for ne02 1:12
% for ne02 13:24
for count_channel = 1:length(limit_channels)

        count_channels = limit_channels(count_channel);
        % --- bin data ---
        count_step = 1;
        
        if if_norm == 1
            
        % default percentile is 10 for df/f normalization
        Fb(count_channel) =  prctile(x(:,count_channels), p);
        
            for o = 1:bin_size:k-bin_size

                y(count_step,count_channel) = (mean(x(o:o+bin_size,count_channels),1)-Fb(count_channel))/Fb(count_channel);
                count_step = count_step+1;
            end
            
        elseif if_norm == 2
            
        Fmin(count_channel) =  prctile(x(:,count_channels), 1);
        Fmax(count_channel) =  prctile(x(:,count_channels), 99);
        
            for o = 1:bin_size:k-bin_size

                y(count_step,count_channel) = (mean(x(o:o+bin_size,count_channels),1)-Fmin(count_channel))/(Fmax(count_channel)-Fmin(count_channel));
                count_step = count_step+1;
            end
            
        
        elseif if_norm == 0
            
            for o = 1:bin_size:k-bin_size

                y(count_step,count_channel) = mean(x(o:o+bin_size,count_channels),1);
                count_step = count_step+1;
            end
        
        end % if_norm
        
end % count_channel
