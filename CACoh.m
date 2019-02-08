function [tvals,pvals] = CACoh(data,fmod,dur,prestim,poststim)

% Cerebro-Acoustic Coherence (CACoh) analysis at fmod frequency (Hz)
% Computes the Magnitude Squared Coherence between an auditory stimulus 
% ?here a click train of duration [dur] and frequency [fmod]? 
% and individual preprocessed IEEG brain responses to the stimulus
% 
% inputs:
% data: a fieldtrip structure of preprocessed and epoched iEEG data with following fields: 
% data.label = Mx1 cell-array containing M channels names
% data.trials = Nx1 cell-array of M-by-P matrices with M = channels and P = time points
% data.time = Nx1 cell-array of 1-by-P matrices P = time points 
% data.fsample = sampling rate of the iEEG data
%
% fmod: a vector of auditory stimulus frequencies (at which coherence with iEEG data is measured)
% 
%%%%%
% Ref.: 
% "The sound of salience: how roughness enhances aversion through neural synchronisation"
% LH Arnal, A Kleinschmidt, L Spinelli, A-L Giraud, P Megevand
% DOI:
% Luc Arnal (luc.arnal@unige.ch)
%%%%%


%% create a stimulus train at fmod frequency, with the same sampling rate as data 
% click duration is 1./data.fsample, train duration is 1
stim = [];
chunk = [1, zeros(1,round(data.fsample.*dur./fmod)-1)];
while numel(stim)<data.fsample*dur
    stim = [stim,chunk];
end
stim = stim(1:data.fsample*dur);

t_prestim = find(data.time{1}>prestim,1,'first');
t_poststim = find(data.time{1}>poststim,1,'first');

%% stim/brain CACoh difference between prestim and poststim windows for each channel

for ch=1:length(data.label)
    ch
    coherence=NaN(length(data.trial),1);
    for tr = 1:length(data.trial)
        % compute coherence in the baseline
        xb = data.trial{tr}(ch,t_prestim:t_prestim+data.fsample*dur-1);
        [cohb, ~] =  mscohere(xb, stim, [], [], [1 fmod],data.fsample);
        % compute coherence in the poststim
        x = data.trial{tr}(ch,t_poststim:t_poststim+data.fsample*dur-1);
        [coha, ~] =  mscohere(x, stim, [], [], [1 fmod],data.fsample);
        coherence(tr) = coha(1,2)-cohb(1,2); % poststim minus prestim coherence
    end
    % store statistics (t and p values)
    [~,pval,~,STATS] = ttest(coherence); 
    tvals(ch) = STATS.tstat;
    pvals(ch) = pval;
end
