function AnalyzeChipChip(tf, clust, chipFile, outFile)

load 'DiffExpGenes' 'psGeneSymbols' 'psNames';
load 'ClusterAnalysis' 'bestKidx' 'K';
load 'TimeLaggedCorr2' 'tlcPV' 'tfNames' 'tfInds';

if clust > K && clust < 1
  error 'invalid cluster';
end

kind = find(clust == bestKidx);
N = length(kind);

tfi = find(strcmp(tfNames, tf));
if length(tfi) == 0
  error 'could not find TF gene name in the list of TFs for which TLC data was compiled';
end

clustGeneNames = cell(1,N);
clustGeneNamesUpper = cell(1,N);
clustGeneTLCPV = zeros(1,N);
clustGeneProbesets = cell(1,N);
for i=1:N
  geneind = kind(i);
  clustGeneNames{i} = psGeneSymbols{geneind};
  clustGeneNamesUpper{i} = upper(clustGeneNames{i});
  clustGeneProbesets{i} = psNames{geneind};
  clustGeneTLCPV(i) = tlcPV(tfi,geneind);
end

[tfs, tfps, tgs, tgps, tgcs, tlccs, otls, mss, pvs]= ...
    textread('TFLagToMotifTargets2.tsv', ...
	     '%s %s %s %s %d %f %d %f %f', 'delimiter', '\t');
clustGeneMotifScores = zeros(1,N);

Q = length(tfs);
for i=1:Q
  if strcmp(tfs{i}, tf)
    tg = tgs{i};
    tgmi = find(strcmp(clustGeneNames, tg));
    if length(tgmi) > 0
      clustGeneMotifScores(tgmi) = mss(i);
    end
  end
end

[ccch, ccss, ccse, ccpv, ccpl, ccpr, ccgh, ...
 ccc5, ccx5l, ccx5r, ...
 ccc3, ccx3l, ccx3r] = textread(chipFile, ...
				'%s %s %s %f %s %s %s %s %s %s %s %s %s', ...
				'headerlines', 1, ...
				'delimiter', '\t');


clustChIPPVs = ones(1,N);
C = length(ccch);
for i=1:C
  close5 = upper(ccc5{i});
  if ~strcmp(close5, 'NULL')
    gi = find(strcmp(clustGeneNamesUpper, close5));
    if length(gi) > 0
%      sprintf('got a chip chip hit for gene: %s', clustGeneNamesUpper{gi})
      clustChIPPVs(gi) = min(clustChIPPVs(gi), ccpv(i));
    end
  end
end

fid=fopen(outFile, 'w+');
for i=1:N
  geneName =clustGeneNames{i};
  cgmScore = '';
  if clustGeneMotifScores(i) > 0
    cgmScore = num2str(clustGeneMotifScores(i));
  end
  chipPV = '';
  if clustChIPPVs(i) < 1
    chipPV = num2str(clustChIPPVs(i));
  end
  geneData = sprintf('%s\t%s\t%s\t%s\t%d\t%s\t%f\t%s\n', ...
		     tf, ...
		     psNames{tfInds(tfi)}, ...
		     geneName, ...
		     clustGeneProbesets{i}, ...
		     clust, ...
		     cgmScore, ...
		     clustGeneTLCPV(i), ...
		     chipPV);
  fprintf(fid, geneData);
end
fclose(fid);
  
