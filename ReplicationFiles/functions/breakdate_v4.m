function [dates, stats]=breakdate_v4(y)
    %  bdate1:  max(tstat) (Chow)
    %  bdate2:  bCUSUM (OB original)
    %  bdate3:  bCUSUM (OB exp)
    n=size(y,1);
    y=y-y(1);
    dy=y(2:n)-y(1:n-1);
    w=exp( (0:n-2)'./(n-1)*3);
    w=w/sqrt(w'*w);
    w=cumsum(w);
    sig=std(dy);
    tstat=zeros(n-4,1);
    bCUSUM=tstat;
    bCUSUM1=bCUSUM;
    uu=tstat;
    for i=3:n-2
        a=y(i-1:n-1)'*dy(i-1:n-1);
        b=y(i-1:n-1)'*y(i-1:n-1);
        %u2=dy(i-1:n-1)-y(i-1:n-1)*a/b;

        % alte Version mit "falschen" Summenindizes
        %a0=sum(dy(i+1:n-1));
        % von mir veränderte Version
        a0=sum(dy(i-1:n-1));

        c=y(1:i-1)'*dy(1:i-1);
        d=y(1:i-1)'*y(1:i-1);
        u1=dy(1:i-1)-y(1:i-1)*c/d;
        u2=dy(i-1:n-1)-y(i-1:n-1)*a/b;
        %uu(i-2)=u1'*u1+u2'*u2;
        uu(i-2)=dy(1:i-1)'*dy(1:i-1)+u2'*u2;
        sig=sqrt(uu(i-2)/n);
        tstat(i-2)=a/sqrt(b)/sig;
        bCUSUM(i-2)=a0/sig/(n-i)^0.5;
        bCUSUM1(i-2)=a0/sig/w(n-i)^0.5;
    end

    %stats = [abs(tstat), abs(bCUSUM), abs(bCUSUM1)];
    stats = [tstat, bCUSUM, bCUSUM1];
    
    % expand stats to original array length
    stats = [nan(2,3); stats; nan(2,3)];
    
    % find maximum and the corresponding index
    [maxval, max_id] =max(stats,[], 1);
    dates = max_id;

    % tstat=[abs(tstat) (1:n-4)'];
    % tstat=sortrows(tstat,1);
    % bCUSUM=[abs(bCUSUM) (1:n-4)'];
    % bCUSUM=sortrows(bCUSUM,1);
    % bCUSUM1=[abs(bCUSUM1) (1:n-4)'];
    % bCUSUM1=sortrows(bCUSUM1,1);
    % bdate1=tstat(end,2)+2;
    % bdate2=bCUSUM(end,2)+2;
    % bdate3=bCUSUM1(end,2)+2;

    %bdate1 = find(abs(tstat)==max(abs(tstat)))+2;
    %bdate2 = find(abs(bCUSUM)==max(abs(bCUSUM)))+2;
    %bdate3 = find(abs(bCUSUM1)==max(abs(bCUSUM1)))+2; 
    
    %stats = [tstat, bCUSUM, bCUSUM1]
end