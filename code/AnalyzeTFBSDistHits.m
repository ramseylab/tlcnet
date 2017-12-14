function AnalyzeTFBSDistHist

[fileList] = textread('GFFFileList.txt','%s');
N = length(fileList);
elements = [];
for i=1:N
  [seq, feat, feat2, coord1, coord2, score, strand, junk, junk2]= ...
      textread(fileList{i}, ...
	       '%s %s %s %d %d %f %s %s %s', ...
	       'headerlines', 1, ...
	       'delimiter', '\t');
  medCoord = 0.5*(coord1+coord2) - 2000;
  elements = [elements; medCoord];
end

median(elements)

set(gca, 'FontSize', 12);
xlabel('Location relative to transcription start site (bp)');
ylabel('Number of motif matches');
hist(elements, 500);

