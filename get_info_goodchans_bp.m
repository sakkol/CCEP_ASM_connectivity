function [bp_good_chans,bp_good_chans_idx,bp_chan_classes] = get_info_goodchans_bp(info,bp_chans)
% input is info prepared by create_elecinfo and output is channels
% excluding outofbrain, artifact patient and artifact block.
% DOES NOT REMOVE SEIZURE OR EPILEPTIC CHANNELS
good_chans_idx = info.channelinfo.outofthebrain~=1&...
                 info.channelinfo.artifact_patient~=1&...
                 info.channelinfo.artifact_block~=1;

known_classes = cell([length(info.channelinfo.Label),1]);
known_classes(good_chans_idx) = {'good'};

% now look at each pair
xx = cellfun(@(x) strsplit(x,'-'),bp_chans,'UniformOutput',0);sp_split=cell([length(xx),2]);for z=1:length(xx),sp_split(z,1:2)=xx{z};end
tobeclassified1=sp_split(:,1);
tobeclassified2=sp_split(:,2);

[bp_chan_classes] = assign_classes(tobeclassified1, tobeclassified2, info.channelinfo.Label, known_classes);

bp_good_chans_idx = find(strcmp(bp_chan_classes,'good'));
bp_good_chans = bp_chans(bp_good_chans_idx);

end