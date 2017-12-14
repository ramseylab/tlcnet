function gen=GetStimuli

[gotStrain, gotStim, gotTimept]=textread('UniqueExptNames.txt', ...
				       '%s %s %s', ...
				       'delimiter', '_');

gen = unique(gotStim);

