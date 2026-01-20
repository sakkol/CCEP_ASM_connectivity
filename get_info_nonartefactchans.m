function [nonartefact_chans,nonartefact_idx] = get_info_nonartefactchans(info)
% input is info prepared by create_elecinfo and output is channels
% excluding only outofbrain and artifact patient and artifact block. 
% Simpler than get_info_goodchans.m.
nonartefact_idx = ...
    info.channelinfo.outofthebrain~=1&...
    info.channelinfo.artifact_patient~=1&...
    info.channelinfo.artifact_block~=1;

nonartefact_chans = info.channelinfo.Label(nonartefact_idx);
end