function [amp_total_pow, amp_env_mean, amp_env_std,...
    amp_skew, amp_kurt] = amplitudeFeatures(x)
%Calculates a bunch of amplitude features
%   x - the input signal
%
% output values
% based on the preprint: NEURAL: quantitative features for newborn EEG using Matlab
%
% amp_total_pow - amplitude total power
% amp_env_mean - mean of the envelope
% amp_env_std - standard deviation of the envelope
% amp_skew - skewness of the signal amplitude
% amp_kurt - kurtosis of the signal amplitude


%---------------------------------------------------------------------
% power of signal (amplitude
%---------------------------------------------------------------------
amp_total_pow = mean( abs(x).^2 );
        
%---------------------------------------------------------------------
% mean or SD of envelope
%---------------------------------------------------------------------
env=abs( hilbert(x) ).^2;       
amp_env_mean = mean(env);
amp_env_std = std(env);
        
%---------------------------------------------------------------------
% skew of amplitude
%---------------------------------------------------------------------
amp_skew = abs(skewness(x));
        
%---------------------------------------------------------------------
% kurtosis of amplitude
%---------------------------------------------------------------------
amp_kurt = kurtosis(x);

end

