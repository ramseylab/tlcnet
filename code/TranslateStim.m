function niceStim = TranslateStim(stim)

if 1==strcmp(stim, 'D4')
  niceStim = 'T091317';
else 
  if 1==strcmp(stim, 'PIC')
    niceStim = 'Poly(I:C)';
  else 
    if 1==strcmp(stim, 'PAM3PIC')
      niceStim = 'Pam3CSK4 / Poly(I:C)';
    else 
      if 1==strcmp(stim, 'R848PIC')
	niceStim = 'R848 / Poly(I:C)';
      else 
	if 1==strcmp(stim, 'D4LPS')
	  niceStim = 'T091317 / LPS';
	else 
	  if 1==strcmp(stim, 'D4PIC')
	    niceStim = 'T091317 / Poly(I:C)';
	  else 
	    if 1==strcmp(stim, 'D4PAM3')
	      niceStim = 'T091317 / Pam3CSK4';
	    else 
	      if 1==strcmp(stim, 'LPSVL')
		niceStim = 'LPS';
	      else 
		if 1==strcmp(stim, 'LPSPAM2')
		  niceStim = 'LPS / Pam2CSK4';
		else 
		  if 1==strcmp(stim, 'PAM2')
		    niceStim = 'Pam2CSK4';
		  else 
		    if 1==strcmp(stim, 'PAM3')
		      niceStim = 'Pam3CSK4';
		    else
		      if 1==strcmp(stim, 'CPG')
			niceStim='CpG';
		      else
			niceStim = stim;
		      end
		    end
		  end
		end
	      end
	    end
	  end
	end
      end
    end
  end
end
