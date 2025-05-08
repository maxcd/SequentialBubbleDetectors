function pval=EoS_test(y,B, detrend, twosided)

if nargin > 4
    error('The function EoS_test takes at most 2 input arguments')
end
if nargin < 4
    twosided=0;
end
if nargin < 3
    detrend=0;
end
if nargin < 2
    B=10;
end

   n=size(y,1);
   dy=y(2:n)-y(1:n-1);
   if detrend
       dy = dy-mean(dy);
   end
   n=n-1;
   nB=n-B;
   Smax=zeros(nB,1);
   for j=1:nB
      dyi=dy(j:j+B-1);
      x=dyi;
      Smax(j)=sum(x)/sqrt(x'*x);  
   end

   if twosided
       ind=abs(Smax(1:nB-B))>abs(Smax(nB));
   else
       ind=Smax(1:nB-B)>Smax(nB);
   end
    
   pval=mean(ind);
    