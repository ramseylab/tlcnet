% CombineAllReplicates.m
%
% Stephen Ramsey, Institute for Systems Biology
% Feb. 2006
%
% For a set of M replicates (over K experiments) and N probesets
% contained in a NxM matrix "logProbeSetIntensities", combine all
% replicates into K columns.  Also, any probeset (row) whose name
% begins with "AFFX-" is eliminated from the data matrix.  This
% function is a wrapper that calls the CombineReplicates function
% once for each of the K experiments, and passes only those columns
% from "logProbeSetIntensities" for the particular replicates for
% the k_i'th experiment.  A cell array is returned as described
% here:
% retDat{1}:  a string array of the K experiment names
% retDat{2}:  a vector (length K) of the number of replicates for
%             each experiment
% retDat{3}:  a NxK matrix of the mean values (in log-base-2)
%             for each experiment and each probeset, after
%             replicate combining.
% retDat{4}:  the standard deviation for the array data,
%             as computed from the mean values and "sigma epsilon"

function retDat=CombineAllReplicates(probeSetIntensities, ...
				     experimentNames)

[numProbeSets,numExperiments]=size(probeSetIntensities);
numParams=numProbeSets+2;

% get the list of unique experiment names
uniqueExperimentNames=unique(experimentNames);

numUniqueExperiments=length(uniqueExperimentNames); 

mu=zeros(numProbeSets,numUniqueExperiments);
experiments={};
numReplicates=zeros(1,numUniqueExperiments);

repInd = cell(1,numUniqueExperiments);
stdDevs = zeros(size(mu));

for condInd=1:numUniqueExperiments 

  % get the name for this experiment/strain/time-point
  experimentName=uniqueExperimentNames{condInd}
  
  % obtain a list of all replicates for this experiment/strain/time-point
  replicateIndices=find(strncmp(experimentNames,experimentName,length(experimentName))==1)
  repInd{condInd}=replicateIndices;
  numReplicates(condInd) = length(replicateIndices);
  experiments{condInd}=experimentName;
  replicateProbeSetIntensities = probeSetIntensities(:,replicateIndices);
  
  mu(:,condInd)=mean(replicateProbeSetIntensities,2);
  if length(repInd) > 1
    stdDevs(:, condInd)=std(replicateProbeSetIntensities,0,2);
  else
    stdDevs(:, condInd)=nan;
  end
end

retDat={};
retDat{1}=experiments;
retDat{2}=numReplicates;
retDat{3}=mu;
retDat{4}=stdDevs;

