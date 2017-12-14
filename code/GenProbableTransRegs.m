function GenProbableTransRegs

load 'DiffExpGenesTFStatus';
load 'DiffExpGenes';

N = length(masterInds);
fid=fopen('ProbableTransRegs.txt','w+');
for i=1:N
  if psIsTR(i)
    fprintf(fid,'%s\n', psGeneSymbols{i});
  end
end
fclose(fid);
  