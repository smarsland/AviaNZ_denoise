%#########################################################################################
%#                                                                                      ##
%# Code of the denoising algorithm presented in "Birdsong Denoising Using Wavelets"     ##
%# published in PLOS One Journal.                                                       ##
%#                                                                                      ##
%# You are free to use, change, or redistribute the code in any way you wish for        ##
%# non-commercial purposes, but please maintain the name of the original author.        ##
%# This code comes with no warranty of any kind.                                        ##
%#                                                                                      ##
%# Nirosha Priyadarshani, Stephen Marsland - November 2015                              ##
%#                                                                                      ##
%#########################################################################################

% Load the birdsongs to be denoised
d = dir(['Primary dataset\kiwi\male\*.wav']);
if isempty(d(:,1)), return, end

for i = 1:length(d) 
    [y,fs] = audioread(['Primary dataset\kiwi\male\' d(i).name]);
    %Find the best decomposition level
    L=BestLevel(y);
    %Generate the wavelet tree using wavelet packet decomposition
    wpt = wpdec(y,L,'dmey');
    %Calculate the detail coefficients on level 1
    det1=wpcoef(wpt,2);
    %Calculate threshold
    sigma=median(abs(det1))/0.6745;
    thr=4.5*sigma;
    %Generate denoised signal
    yd=wpdencmp(y,'s',L,'dmey','threshold',thr,0);
    audiowrite(['Primary dataset\kiwi\male\D\' d(i).name],yd,fs)
    %Band-pass filtering
    D = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',100,500,8000,8400,80,1,80,fs);
    H = design(D,'equiripple'); 
    ydf=filter(H,yd);
    audiowrite(['Primary dataset\kiwi\male\DF\' d(i).name],ydf,fs)
    yf=filter(H,y);
    audiowrite(['Primary dataset\kiwi\male\F\' d(i).name],yf,fs)
end

% Plot the spectrums of last example
windowsize = 128;
window = hanning(windowsize);
nfft = windowsize;
noverlap = windowsize-1;
x=0:(length(y)/fs)/(length(y)-1):length(y)/fs;
clims = [0 75];
[S,F,T] = spectrogram(y,window,noverlap,nfft,fs);
subplot(4,2,1); plot(x',y);axis tight; xlabel('Time (secs)'); title('Original Song');
subplot(4,2,2); imagesc(T,F,log10(abs(S)));
set(gca,'YDir','Normal');xlabel('Time (secs)');ylabel('Freq (Hz)');title('Original Song');

[S,F,T] = spectrogram(yf,window,noverlap,nfft,fs);
subplot(4,2,3); plot(x',yf);axis tight; xlabel('Time (secs)'); title('Band-passed Song');
subplot(4,2,4); imagesc(T,F,log10(abs(S)));
set(gca,'YDir','Normal');xlabel('Time (secs)');ylabel('Freq (Hz)');title('Band-passed Song');

[S,F,T] = spectrogram(yd,window,noverlap,nfft,fs);
subplot(4,2,5); plot(x',yd);axis tight; xlabel('Time (secs)'); title('Denoised Song');
subplot(4,2,6); imagesc(T,F,log10(abs(S)));
set(gca,'YDir','Normal');xlabel('Time (secs)');ylabel('Freq (Hz)');title('Denoised Song');

[S,F,T] = spectrogram(ydf,window,noverlap,nfft,fs);
subplot(4,2,7); plot(x',ydf);axis tight; xlabel('Time (secs)'); title('denoised and band-pass filtered Song');
subplot(4,2,8); imagesc(T,F,log10(abs(S)));
set(gca,'YDir','Normal');xlabel('Time (secs)');ylabel('Freq (Hz)');title('denoised and band-pass filtered Song');

colormap(flipud(gray));