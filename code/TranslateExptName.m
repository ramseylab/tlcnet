function [strain, stim, time] = TranslateExptName(exptName)

if ~iscell(exptName)
  [strain, stim, time]=strread(exptName, '%s_%s_%d');
else
  M = length(exptName);
  strain = cell(1,M);
  stim = cell(1,M);
  time = zeros(1,M);
  for j=1:M
    [jstrain, jstim, time(j)]=strread(exptName{j}, '%s %s %d', ...
				      'delimiter', '_');
    strain{j}=jstrain{1};
    stim{j}=jstim{1};
  end
end

