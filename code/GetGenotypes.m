function gen=GetGenotypes

[gotStrain, gotStim, gotTimept]=textread('UniqueExptNames.txt', ...
				       '%s %s %s', ...
				       'delimiter', '_');

gen = unique(gotStrain);

