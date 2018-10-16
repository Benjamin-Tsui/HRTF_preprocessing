function [frequencyBands] = frequency_bnad_calculation(start_freq, end_freq, bandwidth)
%FREQUENCY_BNAD_CALCULATION Summary of this function goes here
% calculate the cut-off frequency and centre frequency for ovctave or 1/3
% octave band

% start_freq = 20%
% end_freq = 20000%
%bandwidth = 'third-octave';  % or 'octave'

% return a array for frequency bands frequency 
% (1,:) = band center frequency (2,:) = band start frequency (3,:) = band end frequency 


if strcmp(bandwidth,'third-octave') == 1
    start_coeff = round(log2(start_freq/1000)*3);
    end_coeff = round(log2(end_freq/1000)*3);
    fcentre = 10^3 * (2 .^ ((start_coeff:end_coeff)/3));     % 250Hz to 20kHz
    fd = 2^(1/6);
    fupper = fcentre * fd;
    flower = fcentre / fd;
elseif strcmp(bandwidth,'octave') == 1
    start_coeff = round(log2(start_freq/1000));
    end_coeff = round(log2(end_freq/1000));
    fcentre = 10^3 * (2 .^ (start_coeff:end_coeff));        % 250Hz to 16kHz
    fd = 2^(1/2);
    fupper = fcentre * fd;
    flower = fcentre / fd;
else
    error('please enter "octave" or "third-octave" in bandwidth') 
end    
frequencyBands = zeros(length(fcentre), 3);
frequencyBands(:, 1) = flower';
frequencyBands(:, 2) = fcentre';
frequencyBands(:, 3) = fupper';
% cut-off and centre frequency matrix
% column 1 = lower cut-off frequency 
% column 2 = centre frequency
% column 3 = upper cut-off freqency
% row = band number

end

