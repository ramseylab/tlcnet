function GenTimeCourseTable

strains = {'WT';
	   'WT';
	   'WT';
	   'WT';
	   'WT';
	   'ATF3';
	   'ATF3';
	   'ATF3';
	   'ATF3';
	   'CREM';
	   'CREM'};
fullStrains = {'Wild-type';
	       'Wild-type';
	       'Wild-type';
	       'Wild-type';
	       'Wild-type';
	       'Atf3(-/-)';
	       'Atf3(-/-)';
	       'Atf3(-/-)';
	       'Atf3(-/-)';
	       'Crem(-/-)';
	       'Crem(-/-)'};

stims = {'LPS';
	 'PIC';
	 'PAM2';
	 'PAM3';
	 'R848';
	 'CPG';
	 'LPS';
	 'PAM2';
	 'PIC';
	 'LPS';
	 'PIC'};
fullStims = {'LPS';
	     'poly I:C';
	     'Pam2CSK4';
	     'Pam3CSK4';
	     'R848';
	     'CpG';
	     'LPS';
	     'Pam2CSK4';
	     'poly I:C';
	     'LPS';
	     'poly I:C'};
	   
allTimes = cell(0,0);
T = length(stims);
for t=1:T
  [cols, allTimes{t}]=GetColsets(strains{t}, ...
				 stims{t}, ...
				 0);
  allTimes{t}=sort(allTimes{t});
end

fid=fopen('TimeCourseTable.tsv','w+');
for t=1:T
  times = allTimes{t};
  P=length(times);
  timesStr = '';
  for p=1:P
    timesStr = [timesStr sprintf('%d', times(p))];
    if p < P
      timesStr = [timesStr ', '];
    end
  end
 fprintf(fid, '%s\t%s\t%s\n', ...
	    fullStrains{t}, ...
	    fullStims{t}, ...
	    timesStr);
end

fclose(fid);