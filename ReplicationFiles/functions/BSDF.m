function [tstat, adf]  = BSDF(y,p,det,r2)
%--------------------------------------------------------------
% PURPOSE: Performs the backward DF test as in Phillips et al (2015)
%--------------------------------------------------------------
% USAGE: tstat = BSDF(y,r2,d)
% where: y  = a time series vector of dimensions (T,1)
%        r2 = fixed  point of the sample
%        d  = optional scalar indicating whether y has been detrended
%             (d=1 in the case with detrending)
%--------------------------------------------------------------
% RETURNS:
%        tstat = scalar value of the test-statistic
%---------------------------------------------------------------


if nargin > 4
    error('The function supDF takes at most 4 input arguments.')
end

T=length(y);

    
%r0 = floor(0.1*T);
r0 = floor((0.01+1.8/sqrt(T)) *T);
r2 = floor(r2*T);       % fixed endpoint
adf=zeros(r2-r0,1);     % Sequential DF-Test
for r1=1:r2-r0
     adf(r1) = ADFstat(y(r1:r2), p-1, det);
end
tstat=max(adf);

end