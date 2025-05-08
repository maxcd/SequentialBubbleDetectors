function [gdsf, bdsf]  = GSDF(y)
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


if nargin > 1
    error('the function supDF takes at most 2 input arguments')
end

T=length(y);
    
%r0 = floor(0.1*T);
r0 = floor((0.01+1.8/sqrt(T)) *T);

% BSDF for a given r2
bdsf = zeros(T-r0+1,1);
for r2 = r0+1:T
    %r2 = floor(r2*T);
    adf=zeros(r2-r0+1,1);     % Sequential DF-Test  r

    for r1=1:r2-r0
         adf(r1) = ADFstat(y(r1:r2), p-1, det);
    end

    % for r1=1:r2-r0+1
    %      z=y(r1+1:r2);
    %      dz=y(r1+1:r2)-y(r1:r2-1);
    %      b=z\dz;
    %      u=dz-z*b;
    %      sig2=u'*u/(r2-r1-2);  
    %      adf(r1)=b*norm(z)/sqrt(sig2);
    % end;
    bsdf_stat=max(adf);
    bdsf(r2-r0+1) = bsdf_stat;
end
gdsf = max(bdsf);


%fprintf('The test statistic equals ')
end