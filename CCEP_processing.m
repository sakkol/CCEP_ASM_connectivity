function CCEP_processing(Sbj_Metadata, curr_block)
% This is the function that carries out the automatic part of CCEP
% processing.
%
% Input
%    Sbj_Metadata, curr_block - which subject and which block to run
%                   analyses on
% 
% Serdar Akkol, UAB, January 2024.

%% Loading variables and files
% Check if files exist and load them
info_file = fullfile(Sbj_Metadata.iEEG_data,curr_block, [curr_block '_info.mat']);
ecog_var = fullfile(Sbj_Metadata.iEEG_data,curr_block, [curr_block '_ecog_bp.mat']);

% Check if ecog file exists. If not then throw an error
if isfile(ecog_var)
    xx=load(ecog_var);ecog=xx.(['ecog_bp']);clear xx
else
    error('Error! Expected ecog file at %s', ecog_var);
end

% Check if info file exists. If not then throw an error
if isfile(info_file)
    load(info_file,'info')
    events = info.events;
    timeSecs=[];
    for t = 1:length(trials2run)
       timeSecs = [timeSecs, events.indStims{trials2run(t)}];
    end
    timeSecs=timeSecs';
else
    error('Error! Expected info file at %s', info_file);
end

% create the CCEP folder if not created yet
CCEPdir = fullfile(Sbj_Metadata.iEEG_data,curr_block,'CCEPdir');
if ~isfolder(CCEPdir); mkdir(CCEPdir); end

% The name of the file that will contain CCEPs data for analysis
cceps_filename = fullfile(CCEPdir, [curr_block '_ccep_bp_' events.StimCh1{trials2run(1)} '-' events.StimCh2{trials2run(1)} '.mat']);

%% Main Processing
% First resample for saving space and memory
cfg = [];
cfg.resamplefs = 1000;
ecog.ftrip = ft_resampledata(cfg,ecog.ftrip);
fs = ecog.ftrip.fsample;

% define pre and post (beginning and the end of each trial)
prestim_s  = 0.75;
poststim_s = 1.25;

% Make trial structure
trl           = [];
trl(:,1)      = floor( timeSecs*fs - fs*prestim_s );
trl(:,2)      = floor( timeSecs*fs + fs*poststim_s);
trl(:,3)      = floor( -prestim_s*fs );

% Epoch
cfg      = [];
cfg.trl  = trl;
CCEPs    = ft_redefinetrial(cfg,ecog.ftrip);
% replace NaNs with zeros (this happens when stimulation starts in first 3 seconds of recording)
CCEPs.trial{1}(isnan(CCEPs.trial{1})) = 0;
CCEPs.trial{2}(isnan(CCEPs.trial{2})) = 0;
CCEPs.trial{3}(isnan(CCEPs.trial{3})) = 0;

% define preproc parameters
cfg_preproc                     = [];
cfg_preproc.channel             = 'all';
% cfg_preproc.padding             = 4;
% cfg_preproc.demean              = 'yes';
% cfg_preproc.baselinewindow      = [-.5 -.05];
% cfg_preproc.detrend             = 'yes';
cfg_preproc.dftfilter           = 'yes';
cfg_preproc.dftfreq             = [60 120 180];
CCEPs                           = ft_preprocessing(cfg_preproc, CCEPs);

% using ft_timelockbaseline for baseline correction
cfg              = [];
cfg.baseline     = [-.215 -.015];
cfg.channel      = 'all';
cfg.parameter    = 'trial';
CCEPs            = ft_timelockbaseline(cfg, CCEPs);

% compute z-score but keep raw values ????????????????????
%         data.trial_raw = data.trial;
%         data.trial = cellfun(@(x) zscore(x,1,2),data.trial,'UniformOutput',0);

% compute ccep trials
cfg             = [];
cfg.keeptrials  = 'yes';
CCEPs           = ft_timelockanalysis(cfg,CCEPs);

%% Save
fprintf('\nSaving to %s\n',cceps_filename)
save(cceps_filename,'CCEPs','-v7.3');

end