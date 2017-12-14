function tgtStr = TranslateStrain(gtStr)

if ~strcmp(gtStr, 'WT') && ~strcmp(gtStr, 'TRIF')
  tgtStr = ['{\it ' gtStr(1) lower(gtStr(2:length(gtStr))) '}^{(-/-)}'];
else
  if strcmp(gtStr, 'TRIF')
    tgtStr = '{\it Ticam1}^{(Lps2/Lps2)}';
  else
    tgtStr = 'WT';
  end
end
