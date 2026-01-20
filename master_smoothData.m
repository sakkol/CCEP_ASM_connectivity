function [smoothed_ftrip] = master_smoothData(ftrip,smooth_factor)
% This function runs smoothing using gaussian convolution on time-series or
% on time-frequency data. If TF data is given, freq_OI is needed (default
% is 70-170Hz HFA). smooth_factor is to create the window size for gaussian
% kernel (default is 0.1).
%
% Serdar Akkol, Human Brain Mapping Lab
% May 2022

%% Check inputs
if ~isfield(ftrip, {'powspctrm','trial'})
    error('Input ftrip is not recognized. Make sure it has powspctrm or trial field.')
elseif isfield(ftrip, 'powspctrm')
    to_conv = ftrip.powspctrm;
    fprintf('Input is assumed to be 4D time-frequency data (trial x channel x frequency x time).\n')
elseif isfield(ftrip, 'trial')
    to_conv = ftrip.trial;
    fprintf('Input is assumed to be 3D time-series data (trial x channel x time).\n')
end
    

% smooth_factor
if ~exist('smooth_factor','var') || isempty(smooth_factor)
    smooth_factor = 0.1;
end
winSize = floor((1/(ftrip.time(2)-ftrip.time(1)))*smooth_factor);
gusWin= gausswin(winSize)/sum(gausswin(winSize));

%% Do the job based on size of input
if length(size(to_conv)) == 2
    error('This input is not recognized, what is it?')
elseif length(size(to_conv)) == 3
    conved = zeros(size(to_conv));
    for rpt = 1:size(to_conv,1)
        for chan = 1:size(to_conv,2)
            conved(rpt,chan,:) = convn(squeeze(to_conv(rpt,chan,:)),gusWin,'same');
        end
    end
    
elseif length(size(to_conv)) == 4
    conved = zeros(size(to_conv));
    for rpt = 1:size(to_conv,1)
        for chan = 1:size(to_conv,2)
            for freq = 1:size(to_conv,3)
                conved(rpt,chan,freq,:) = convn(squeeze(to_conv(rpt,chan,freq,:)),gusWin,'same');
            end
        end
    end
    
end

%% Get the output ready
smoothed_ftrip = ftrip;
if isfield(ftrip, 'powspctrm')
    smoothed_ftrip.powspctrm = conved;
elseif isfield(ftrip, 'trial')
    smoothed_ftrip.trial = conved;
end

end