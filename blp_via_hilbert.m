function BLP=blp_via_hilbert(data,freqs,Fs)
% function BLP=blp_via_hilbert(data,freqs,Fs)
%
% compute bandpass power of a specific electrode
%
% data is a row vector of data
% freqs is a vector defining band boundaries (e.g 50:10:150)
% Fs is sampling rate


%means_sum=0;
if size(data,1)>1
   error('data needs to be a row vector'); 
end

n_band=length(freqs)-1;
n_tpt=size(data,2);
BLP=zeros(n_band,n_tpt);

for j=1:n_band,
    lowF=freqs(j);
    highF=freqs(j+1);
    
    if lowF~=0 && highF~=0 % a workaround to avoid the 'singualr matrix' error                                    
            [B A]=butter(4,[lowF highF]*2/Fs);
            BLP(j,:)=abs(hilbert(filtfilt(B,A,data)));
    else
        BLP(j,:)=abs(hilbert(eegfilt(data,Fs, lowF, highF)));
    end
end

BLP=bsxfun(@rdivide,BLP,mean(BLP,2)); %this averages across time points but within freq bands
%Thus you normalize the time series in each band by dividing by its mean

BLP=mean(BLP,1).*100; %percent of mean signal, this averages across freq bands