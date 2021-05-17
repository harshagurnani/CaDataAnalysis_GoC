function Hd = whisk_setpoint_FIR( Fs )
%GETFILTER Returns a discrete-time filter object.


Fpass = 1;  % Passband Frequency
Fstop = 5;   % Stopband Frequency
Apass = 1;     % Passband Ripple (dB)
Astop = 60;    % Stopband Attenuation (dB)


% Fpass = 0.05;  % Passband Frequency
% Fstop = 0.5;   % Stopband Frequency
% Apass = 1;     % Passband Ripple (dB)
% Astop = 60;    % Stopband Attenuation (dB)

if isempty(Fs),
Fs    = 300;   % Sampling Frequency
end

h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, Fs);

Hd = design(h, 'equiripple','MinOrder', 'any', 'StopbandShape', 'flat');


end