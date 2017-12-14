function GenUniqueExptNamesFile

load 'DiffExpGenes' 'uniqueExptNames' 'numReplicates';

M=length(uniqueExptNames);

fid=fopen('UniqueExptNames.tsv', 'w+');
for i=1:M
  fprintf(fid, '%d\t%s\t%d\n', i, ...
	  uniqueExptNames{i}, ...
	  numReplicates(i));
end
fclose(fid);


