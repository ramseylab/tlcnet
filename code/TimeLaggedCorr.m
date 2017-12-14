function tlc = TimeLaggedCorr(regDat, targDat, shifts, cols, times)
% Computes the time-lagged correlation (TLC) coefficient for a large
% number of pairs of transcription factors and potential target genes.
% The entry point for the code is the function TimeLaggedCorr.m. The
% first two arguments to this function are the matrices "regDat" and
% "targDat".  The "regDat" matrix contains the gene expression
% measurements of the possible transcriptional regulators, one row per
% gene and one column per measurement (multiple time-courses can be
% contained in this matrix, and the time points & conditions do not
% need to be consecutive).  The "targDat" matrix contains the gene
% expression measurements of the possible target genes, with one gene
% per row and one condition/measurement per column. The column
% ordering of the "targDat" matrix must be consistent with the column
% ordering of the "regDat" matrix.  The "shifts" argument is a vector
% that specifies the various possible time lags for which you wish to
% compute the TLC coefficient. The "cols" argument is a cell array,
% with each element of the cell array corresponding to a specific
% time-course (each time-course is represented by one element in the
% cell array).  For each time-course study, the corresponding element
% in the "cols" cell array contains an array of indices that specify
% the columns of the "regDat" and "targDat" matrices that correspond
% to the measurements in that time-course (in chronological order).
% The "times" argument is a cell array that has the same length and
% ordering as the "cols" argument, and each element of this cell array
% contains an array that specifies the time points for each
% measurement in each time course study (in the same units as the time
% values in the "shifts" vector).  The function returns a
% three-dimensional matrix (of the number of elements of "shifts" x
% the number of potential regulators x the number of potential target
% genes) containing the time-lagged correlation coefficients.

NR = size(regDat, 1);
NT = size(targDat, 1);
M = size(regDat, 2);

if M ~= size(targDat, 2)
  error 'targDat matrix should have the same number of columns as the regDat matrix';
end

T = length(shifts);

allCols = [];
K = length(cols);
if K ~= length(times)
  error 'times and cols should be the same length';
end

for k=1:K
  ck=cols{k};
  ckl=length(ck);
  allCols = [allCols ck(1:(ckl-1))'];
end

tlc = zeros(T, NR, NT);

for t=1:T
  shiftTime = shifts(t);
  
  for i=1:NT
    tlcvec = TimeLaggedCorrList(shiftTime, ...
				regDat, ...
				targDat(i,:), ...
				cols, ...
				times);
    tlc(t, :, i) = tlcvec(:);
  end
end

