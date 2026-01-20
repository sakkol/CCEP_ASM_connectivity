function [ccep_rms] = ccep_rms(avg_amp,time_line,part)

if ~exist('part','var') || isempty(part)
    part = [];
end
if isempty(part) || strcmp(part,'all')
    ccep_rms = rms(avg_amp(time_line>0.010 & time_line<0.300));
elseif strcmp(part,'N1')
    ccep_rms = rms(avg_amp(time_line>0.010 & time_line<0.050));
elseif strcmp(part,'N2')
    ccep_rms = rms(avg_amp(time_line>0.050 & time_line<0.300));
end

end