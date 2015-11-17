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

function [L]=BestLevel(y)
E=wentropy(y,'shannon');
prefered=1;     % prefer to further decompose
L=1;            % initially level is 1
while ((prefered==1)|| (14<L))
    wpt = wpdec(y,L,'dmey');    % decompose (level L)
    % calculate the entropy E of each leaf node
    last=2^L-1;
    e=zeros(1,last);
    for i=0:last
         N = depo2ind(2,[L i]);
         det = wpcoef(wpt,N);
         e(i+1) = wentropy(det,'shannon');
    end
    if (max(e)> E)
        prefered=0;
        L=L-1;
        return
    end
    E=max(e);
    clear e;
    L=L+1;
end