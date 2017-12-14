function GenTableExperiments

[cols, exptNames, numReps]=textread('UniqueExptNames.tsv','%d %s %d');

strains = {'ATF3';
	   'CREM';
	   'MYD88';
	   'TRIF';
	   'WT'};
strainsFull = {'Atf3(-/-)';
	       'Crem(-/-)';
	       'Myd88(-/-)';
	       'Ticam1(Lps2/Lps2)';
	       'Wild-type'};

stims = {'CPG';
	 'LPS';
	 'LPSPAM2';
	 'PAM2';
	 'PAM3';
	 'PAM3PIC';
	 'PIC';
	 'R848'};
stimsFull = {'CpG';
	     'LPS';
	     'LPS / Pam2CSK4';
	     'Pam2CSK4';
	     'Pam3CSK4';
	     'Pam3CSK4 / poly I:C';
	     'poly I:C';
	     'R848'};

S = length(strains);
A = length(stims);

fid=fopen('ExperimentsTable.tsv','w+');
for s=1:S
  strain = strains{s};
  [cols, times]=GetColsets(strain, 'UNSTIM', 0);
  fprintf(fid, '%s\t%s\t%d\t%d\n', ...
	  strainsFull{s}, ...
	  'unstimulated', ...
	  0, ...
	  numReps(cols(1)));
  for a=1:A
    stim = stims{a};
    [cols, times]=GetColsets(strain, ...
			     stim, ...
			     1);
    T = length(times);
    for t=1:T
      fprintf(fid, '%s\t%s\t%d\t%d\n', ...
	      strainsFull{s}, ...
	      stimsFull{a}, ...
	      times(t), ...
	      numReps(cols(t)));      
    end
  end
end
fclose(fid);