function [rEEG_CV, rEEG_asymmetry] = rangeEEG(x,Fs)
%-------------------------------------------------------------------------------
% rEEG: range EEG as defined in [1]
% [1] D O’Reilly, MA Navakatikyan, M Filip, D Greene, & LJ Van Marter (2012). Peak-to-peak
% amplitude in neonatal brain monitoring of premature infants. Clinical Neurophysiology,
% 123(11), 2139–53.
%


% generate rEEG
% 4 seconds and 50% overlap
reeg=gen_rEEG(x,Fs,4,50,'hamm');


%---------------------------------------------------------------------
% coefficient of variation
%---------------------------------------------------------------------
rEEG_CV = std(reeg)/mean(reeg);

%---------------------------------------------------------------------
% measure of assymetry
%---------------------------------------------------------------------
A=median(reeg) - prctile(reeg,5);
B=prctile(reeg,95) - median(reeg);
rEEG_asymmetry =(B-A)/(A+B);

end

function reeg=gen_rEEG(x,Fs,win_length,win_overlap,win_type)
%---------------------------------------------------------------------
% generate the peak-to-peak measure (rEEG)
%---------------------------------------------------------------------
[L_hop,L_epoch,win_epoch]=gen_epoch_window(win_overlap,win_length,win_type,Fs);


N=length(x);
N_epochs=floor( (N-(L_epoch-L_hop))/L_hop );
if(N_epochs<1) N_epochs=1; end
nw=0:L_epoch-1;

%---------------------------------------------------------------------
% generate short-time FT on all data:
%---------------------------------------------------------------------
reeg=NaN(1,N_epochs);
for k=1:N_epochs
    nf=mod(nw+(k-1)*L_hop,N);
    x_epoch=x(nf+1).*win_epoch(:)';

    reeg(k)=max(x_epoch)-min(x_epoch);
end
% no need to resample (as per [1])


% log--linear scale:
ihigh=find(reeg>50);
if(~isempty(ihigh))
    reeg(ihigh)=50.*log(reeg(ihigh))./log(50);
end
end
