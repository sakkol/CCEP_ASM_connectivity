function [cell_splitted] = strsplit_SA(to_split,delimeter)
% To split the bipolar labels into cell, to decrease clutter.

% choose delimeter
if ~exist('delimeter','var') || isempty(delimeter)
    delimeter = '-';
end

xx = cellfun(@(x) strsplit(x,delimeter),to_split,'UniformOutput',0);
cell_splitted=cell([length(xx),2]);
for z=1:length(xx)
    cell_splitted(z,1:2)=xx{z};
end