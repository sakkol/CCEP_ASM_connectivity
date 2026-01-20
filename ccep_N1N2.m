function [peak_amps,peak_times,N1_amp,N1_time,N2_amp,N2_time] = ccep_N1N2(avg_data,time_line)
% Tidying up the CCEP_com_offon to calculate N1 and N2 responses with
% findpeaks.m. N1 is within 10-30ms and N2 is within 80-250ms.

[~,tt] = findpeaks(abs(avg_data));
peak_amps = avg_data(tt);
peak_times = time_line(tt);

% N1:
N1_idx = find(peak_times>0.01 & peak_times<0.03);
all_N1_amps = peak_amps(N1_idx);
[~,xx] = max(abs(all_N1_amps));
N1_amp  = all_N1_amps(xx);
N1_time = peak_times(N1_idx(xx));
if isempty(N1_amp)
    N1_amp=NaN;
    N1_time=NaN;
end

% N2:
N2_idx = find(peak_times>0.085 & peak_times<0.25);
all_N2_amps = peak_amps(N2_idx);
[~,xx] = max(abs(all_N2_amps));
N2_amp  = all_N2_amps(xx);
N2_time = peak_times(N2_idx(xx));
if isempty(N2_amp)
    N2_amp=NaN;
    N2_time=NaN;
end

end