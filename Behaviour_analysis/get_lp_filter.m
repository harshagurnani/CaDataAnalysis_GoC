function [sosbp, gbp] = get_lp_filter( Fs, Wp, Ws )

%     Fs = 300;%length(Whiskers_time_0)/((Whiskers_time_0(end)-Whiskers_time_0(1))/1000);
    Fn = Fs/2;                                              % Nyquist Frequency
    Wp = Wp/Fn;%lowpass: 2%3/Fn;        slow: 0.3/Fn;                               % Passband Frequencies (Normalized)
    Ws = Ws/Fn;%8/Fn;                        1/Fn;                % Stopband Frequencies (Normalized)
    Rp = 10;                                                % Passband Ripple (dB)
    Rs = 50;                                                % Stopband Ripple (dB)
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                         % Filter Order
    [c,b,a] = cheby2(n,Rs,Ws);                              % Filter Design
    [sosbp,gbp] = zp2sos(c,b,a);                            % Convert To Second-Order-Section For Stability
    % WSP = filtfilt(sosbp, gbp, Whiskers_angle_0);
end