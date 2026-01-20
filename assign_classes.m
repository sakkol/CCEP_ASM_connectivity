function [assigned_classes] = assign_classes(tobeclassified1, tobeclassified2, known_names, known_classes, importance)
% This function makes is easier to assign 1 classification to a pair of
% names. For example, this can be used to assign seizure onset zone vs
% epileptogenic zone vs healthy etc to a pair of electrodes depending on
% the known classifications of the pair. 
% So:
% - tobeclassified1 and tobeclassified2 are the pair of electrodes,
% - known_names is the names of each electrode,
% - known_classes is the corresponding classifications for "known_names"
% 
% Has the functionality to assign the value based on what is more
% important. For example, if one of pair is SOZ, other is healthy, no
% matter the orientation of pair (elec1-elec2 vs elec2-elec1), it will
% assign SOZ if "importance" is given. Otherwise, will assign what is the
% first of the electrodes.
%
% Serdar Akkol, MD PhD
% Neurology, UAB
% January, 2024
%
%% Check input
% check length of pairs
if length(tobeclassified1) ~= length(tobeclassified2)
    error('Please check the pairs, they seem to be in different lengths.')
end

% check if each pair is in known_names
if any(~ismember(tobeclassified1,known_names))
    error('%s is not in known_names.',strjoin(tobeclassified1(~ismember(tobeclassified1,known_names)),','))
end
if any(~ismember(tobeclassified2,known_names))
    error('%s is not in known_names.',strjoin(tobeclassified2(~ismember(tobeclassified2,known_names)),','))
end

%% startup
assigned_classes = cell([length(tobeclassified1),1]);

%% Loop the to-be-classifieds

for t = 1:length(tobeclassified1)
    % find what the pairs are assigned
    tmp_class1 = known_classes(strcmp(known_names, tobeclassified1{t}));
    tmp_class2 = known_classes(strcmp(known_names, tobeclassified2{t}));

    % check if both pairs are assigned different things
    if ~strcmp(tmp_class1,tmp_class2) && ~exist('importance','var')
        warning('%s-%s pair is assigned two different classes, assigning the classification of first of the pair as the pair''s class!',tobeclassified1{t},tobeclassified2{t})
        assigned_classes(t) = tmp_class1;
    elseif ~strcmp(tmp_class1,tmp_class2) && exist('importance','var')
        tmp_class1_imp = find(ismember(importance,tmp_class1));
        tmp_class2_imp = find(ismember(importance,tmp_class2));
        if tmp_class1_imp < tmp_class2_imp
            assigned_classes(t) = tmp_class1;
        elseif tmp_class1_imp > tmp_class2_imp
            assigned_classes(t) = tmp_class2;
        end
    elseif strcmp(tmp_class1,tmp_class2)
        assigned_classes(t) = tmp_class1;
    end
    
    % collect output

end

end