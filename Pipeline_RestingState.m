%% Preprocessing the BHA from the average data and correlating between each 
% channel to obtain functional connectivity matrices

data_root = fullfile(cd,'neural_data');
project_name = 'CCEP_PrePost';

for subj = 1:2
    if subj==1, sbj_ID = 'P2';
    elseif subj==2, sbj_ID = 'P4';
    end
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID);

    %% Calculate HFB and slow fluctuations; save channels individually
    load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']))

    freq_range = [1 200];
    nf = length(freq_range(1):3:freq_range(2));
    all_freqs = logspace(log10(freq_range(1)),log10(freq_range(2)),nf);
    foi_hfb = all_freqs(all_freqs>69);

    save_dir = fullfile(Sbj_Metadata.iEEG_data, curr_block,'HFB_files');
    if ~exist(save_dir,'dir'),mkdir(save_dir),end

    for c = 1:length(ecog_avg.ftrip.label)

        fprintf(['\n\nNow, electrode number ' num2str(c) '\n\n'])
        cfg                   = [];
        cfg.channel           = ecog_avg.ftrip.label{c};
        % Wavelet
        cfg.method            = 'wavelet';
        % nf                    = length(freq(1):3:freq(2));
        cfg.foi               = foi_hfb;
        nb_cycles = zeros(size(cfg.foi));
        nb_cycles(cfg.foi<2)                = 4;
        nb_cycles(cfg.foi<4 & cfg.foi>=2)   = 4;
        nb_cycles(cfg.foi<9 & cfg.foi>=4)   = 5;
        nb_cycles(cfg.foi<13 & cfg.foi>=9)  = 6;
        nb_cycles(nb_cycles==0)             = 7;
        cfg.width             = nb_cycles;
        cfg.toi               = 2:0.01:[floor(ecog_avg.ftrip.time{1}(end))-2];
        cfg.keeptrials        = 'yes';
        cfg.keeptaper         = 'yes';
        cfg.output            = 'fourier';
        cfg.pad               = 'nextpow2';
        hfb_ftrip_fourier     = ft_freqanalysis(cfg,ecog_avg.ftrip);

        hfb_ftrip_pow = hfb_ftrip_fourier;
        hfb_ftrip_pow.powspctrm = abs(hfb_ftrip_fourier.fourierspctrm).^2;
        hfb_ftrip_pow = rmfield(hfb_ftrip_pow,'fourierspctrm');

        % Compute baseline statistics (geometric mean for log transform)
        baseline_geomean = exp(mean(log(hfb_ftrip_pow.powspctrm), 2));

        % Log transform normalization for each frequency
        normalized_power = log(hfb_ftrip_pow.powspctrm ./ baseline_geomean);

        % Average across frequencies to get HFB
        hfb = squeeze(mean(normalized_power, 3));
        hfb_ftrip_pow.powspctrm = hfb';

        % Create a new data structure for filtering
        cfg = [];
        cfg.time = {hfb_ftrip_pow.time};
        cfg.trial = {hfb'};
        cfg.fsample = 100;
        cfg.label = hfb_ftrip_pow.label;
        temp_data = ft_preprocessing([], cfg);

        % Bandpass filter to get slow fluctuations
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [0.1 1];
        cfg.bpfiltord = 3;
        cfg.bpfilttype = 'but';  % Use BUT filter for better phase preservation
        cfg.usefftfilt = 'yes';   % Use FFT-based filtering for efficiency
        slow_hfb = ft_preprocessing(cfg, temp_data);

        % Save data to iEEG data folder
        s = [];
        s.hfb_ftrip_pow=hfb_ftrip_pow;
        s.slow_hfb=slow_hfb;
        s.hfb_ftrip_fourier=hfb_ftrip_fourier;
        save(fullfile(save_dir,[ecog_avg.ftrip.label{c} '_HFB.mat']),'hfb_ftrip_pow','slow_hfb','hfb_ftrip_fourier','-v7.3')
        % parsave(fullfile(save_dir,[ecog_avg.ftrip.label{c} '_HFB.mat']),'hfb_ftrip_pow','slow_hfb','hfb_ftrip_fourier')

        clear slow_hfb temp_data
    end

    %% Correlate the slow fluctuations
    load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
    save_dir = fullfile(Sbj_Metadata.iEEG_data, curr_block,'HFB_files');

    ch_count = length(info.channelinfo.Label);
    corrmatrix = zeros([ch_count,ch_count]);
    corrmatrix_pval = zeros([ch_count,ch_count]);

    for c1 = 1:ch_count
        c1_sHFB = load(fullfile(save_dir,[info.channelinfo.Label{c1} '_HFB.mat']),'slow_hfb');

        parfor c2 = 1:ch_count
            disp(['c1: ' num2str(c1) '; c2: ' num2str(c2)])
            c2_sHFB = load(fullfile(save_dir,[info.channelinfo.Label{c2} '_HFB.mat']),'slow_hfb');
            [corrmatrix(c1,c2), corrmatrix_pval(c1,c2)] = corr(squeeze(c1_sHFB.slow_hfb.trial{1})', squeeze(c2_sHFB.slow_hfb.trial{1})');
        end

    end

    save(fullfile(save_dir,'corr_matrix.mat'),'corrmatrix','corrmatrix_pval')

end
