% arrayDat:  an NxQ matrix of expression data, where N is the
% number of genes to be tested, and Q is the number of experiments.
% arrayNames:   not used; just pass em
%
%


function [diffExpGenes, ...
	  pValues, ...
	  pValueCutoffs, ...
	  tcExptNames]=...
         GetDiffExpGenes(arrayDat, ...
			 arrayNames, ...
			 exptNames, ...
			 stims, ...
			 strains, ...
			 numIter, ...
			 sl, ...
			 numBins, ...
			 fdr)


% numIter:  number of bootstrap iterations
[N, Q]=size(arrayDat);

timeCourseExperimentsNames = {};
timeCourseExperimentsTimes = {};

numStims = length(stims);
numStrains = length(strains);

tcExptNames = {};

exptCtr = 0;
for i = 1:numStrains
  strainName = strains{i};
  for j=1:numStims
    stimName = stims{j};
    
    [cols, exptTimes]=GetColsets(strainName, stimName, 0);
    
    % get the time points available for this combination of strain
    % and stimulus
    exptTimes = sort(exptTimes);
    
    T = length(exptTimes);
    
    timeCourseExptNames = {};
    timeCourseExptNames{1} = sprintf('%s_UNSTIM_0000', strainName);
    
    for k=2:T
      timeCourseExptNames{k} = sprintf('%s_%s_%0.4d', strainName, stimName, ...
				       exptTimes(k));
    end

    exptCtr = exptCtr + 1;
    
    tcExptNames{exptCtr} = sprintf('%s_%s', strainName, stimName);
    timeCourseExperimentsNames{exptCtr} = timeCourseExptNames;
    timeCourseExperimentsTimes{exptCtr} = exptTimes';
  end
end

M = length(timeCourseExperimentsNames);

pValues = ones(N,M);

tic

% go through each time-course experiment, one at a time
for j=1:M
  sprintf('working on experiment %d out of %d', j, M)
  
  % get the name of this particular time-course experiment
  timeCourseExptNames = timeCourseExperimentsNames{j};
  
  % get the time points associated with this time-course experiment
  timeCourseExptTimes = timeCourseExperimentsTimes{j};
  
  % get the number of time points associated with this time-course experiment
  T = length(timeCourseExptNames);
  
  % build vectors to contain the indices and time points associated
  % with individual arrays for this experiment
  timeCourseExptArrayInds = [];
  timeCourseExptArrayTimes = [];
  
  if T < 4
    timeCourseExptNames
    error 'insufficient number of samples for experiment';
  end
  
  % choose the order of the polynomial to use
  P = min(T - 2, 4);

  % build a complete vector of array samples for this time-course experiment
  for k=1:T
    exptName = timeCourseExptNames{k};

    exptInds = find(strcmp(exptName, exptNames))';
    
    timeCourseExptArrayInds = [timeCourseExptArrayInds exptInds];
    
    exptTime = timeCourseExptTimes(k);
    timeCourseExptArrayTimes = [timeCourseExptArrayTimes exptTime*ones(size(exptInds))];
  end
  
  S = length(timeCourseExptArrayInds);
  
  % get the expression measurements for this
  exptDat = arrayDat(:,timeCourseExptArrayInds);
  
  timeCourseExptArrayTimesMean = mean(timeCourseExptArrayTimes);
  timeCourseExptArrayTimesStd = std(timeCourseExptArrayTimes);
  timeCourseExptArrayTimesTransf = (timeCourseExptArrayTimes - ...
				    timeCourseExptArrayTimesMean)/ ...
      timeCourseExptArrayTimesStd;
  
  for i=1:N
    exptDatProbe = exptDat(i,:);
    
    % perform polynomial fit for alternative hypothesis
    [pfijalt,Salt] = polyfit(timeCourseExptArrayTimesTransf, exptDatProbe, P);

    % get norm of residuals for alternative hypothesis
    nalt = Salt.normr;

    % get norm of residuals for null hypothesis
    nnull = std(exptDatProbe)*sqrt(S-1);
    
    % get residuals from the alternative model
    altRes = polyval(pfijalt, timeCourseExptArrayTimesTransf)- ...
	     exptDatProbe;

    meanExptDatProbe = mean(exptDatProbe);
    
    % get the residuals from the null model
    nullRes = exptDatProbe - meanExptDatProbe;
    
    if S ~= length(altRes)
      error 'incorrect vector size';
    end
    
    sampleInds = ceil(S*rand(numIter, S));
    
    % sample the residuals of the alternative model, "numIter" times
    nnullsamp = repmat(nullRes, numIter, 1) + altRes(sampleInds);

    % create a matrix of the null model fit vals
    nullmodelfitvalsmat = repmat(meanExptDatProbe, numIter, S);
    
    % add the sampled residuals to the null model fit
    nullmodelvalssamp = nullmodelfitvalsmat + nnullsamp;
    
    nullsamp = sqrt(S-1)*std(nullmodelvalssamp,[],2);
    
    Fsamp = (repmat(nnull,numIter,1)-nullsamp)./nullsamp;
    Fi = (nnull-nalt)/nalt;
    
    pv = GaussKernDensCDFTwoTailed(Fsamp', Fi, sl, numBins);
    pValues(i,j)=pv;
  end
end

diffExpGenes = zeros(1,N);

pValueCutoffs = zeros(1,M);
for j=1:M
  pvc = EstimateFDR(pValues(:,j)', fdr);
  pValueCutoffs(j)=pvc;
  ind = find(pValues(:,j) <= pvc);
  diffExpGenes(ind) = 1;
end

elapsedTime = toc;

timePerGene = elapsedTime / (M*N);

sprintf('elapsed time per gene, per experiment: %f sec', timePerGene)


