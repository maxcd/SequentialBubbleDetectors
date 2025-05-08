function [bdate1,bdate2,bdate3, stats]=breakdate_v3b(y, n0)
    %  bdate1:  max(tstat)/Chow
    %  bdate2:  SSR/max(logLikelihood)
    %  bdate3:  bCUSUM
    n=size(y,1);
    %y=y-mean(y);
    y=y-y(1);
    dy=y(2:n)-y(1:n-1);
    %sig=std(dy)

    %tstat=zeros(n-4,1);
    tstat=nan(n,1);
    bCUSUM=tstat;
    uu=tstat;
    uu_R = uu;
    sigma = zeros(n,1);
    alist = zeros(n, 1);
    
    nmin = 2; % miminum window length
    if nargin<2
        n0 = nmin; 
    end
    
    for i=n0+1:n-nmin

        % OLS in the second part of the sample
        y2  = y(i-1:n-1);
        dy2 = dy(i-1:n-1);
        XX2 = y2'*y2; 
        XY2 = y2'*dy2;
        %b2 = nom2/denom2;
        b2  =  y2\dy2;%XY2/XX2;
        %a=y(i-1:n-1)'*dy(i-1:n-1);
        %b=y(i-1:n-1)'*y(i-1:n-1);
        a = XY2; b = XX2;
        
        % Data in the first part of the sample
        y1  = y(1:i-2);
        dy1 = dy(1:i-2);
        b1  = y1\dy1;

        %u2=dy(i-1:n-1)-y(i-1:n-1)*a/b;
        %a0=sum(dy(i+1:n-1));
        a0=sum(dy(i-1:n-1));
        c=y(1:i-1)'*dy(1:i-1);
        d=y(1:i-1)'*y(1:i-1);
        
        % Unrestricted resiuals in second part of the sample
        %u2=dy(i-1:n-1)-y(i-1:n-1)*a/b;
        u2 = dy2 - y2*b2;

        % Restricted sigma and RSS: Imposing the H0 in the first part of the sample
        u1_R=dy(1:i-2);%-y(1:i-1)*c/d;
        uu_R(i)=u1_R'*u1_R+u2'*u2;
        sig_R=sqrt(uu_R(i)/n);
        
        % unrestricted sigma and RSS
        u1=dy1-y1*b1;%c/d;
        uu(i)=u1'*u1+u2'*u2;
        sig=sqrt(uu(i)/n);

        % compute test statistics
        tstat(i)=a/sqrt(b);
        bCUSUM(i)=a0/sqrt(n-i+1);
        sigma(i) = sig_R;
        alist(i) = a0;
    end

    uu_ML = uu_R;
    % tstat=[abs(tstat) (1:n-4)'];
    % tstat=sortrows(tstat,1);
    % uu=[uu (1:n-4)'];
    % uu=sortrows(uu,1);
    % bCUSUM=[abs(bCUSUM) (1:n-4)'];
    % bCUSUM=sortrows(bCUSUM,1);
    % 
    % bdate1=tstat(end,2);
    % bdate2=uu(1,2);
    % bdate3=bCUSUM(end,2);

    bdate1 = find(tstat==max(tstat));%+n0;
    %bdate1 = find(tstat==max(tstat))+2;
    bdate2 = find(uu_ML==min(uu_ML));%+n0;
    bdate3 = find(bCUSUM==max(bCUSUM));%+n0; 

    stats = [tstat, -uu_ML, bCUSUM];
    %stats = [nan(n0,3); stats; nan(2,3)];
    %[maxval, max_id] =max(stats,[], 1);

    
end