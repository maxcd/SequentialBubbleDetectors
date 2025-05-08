function y = simSingleBubble(T, rho, sig2, lam0, lam1, y0, mu, beta, B0, allow_neg_bubbles, eps)

    
    if nargin < 10
         %eps = randn(T,1);
         allow_neg_bubbles = 0;
    end

     if nargin < 11
         eps = randn(T,1);
         %allow_neg_bubbles = 0;
     end

     if nargin < 9
         B0=[];
     end

    t0 = floor(T*lam0);

    t1 = floor(lam1*T);
    n0 = t0;            % Random Walk perioden
    n1 = t1-t0;         % explosive perioden
    if lam0==0
        t0=1;
        n0=0;
        n1=T;
    end
    
    rhot = [ones(max(t0-1, 0), 1); repmat(rho,min(t1-t0, T-1), 1)];
    H = speye(t1) - sparse(2:t1, 1:t1-1, rhot, t1, t1);

    HB = speye(n1) - sparse(2:n1, 1:n1-1, rho, n1, n1);


       pos=0;
       while (pos==0)
            det = ones(T,1)*mu + (1:T)'*beta;
            u = sqrt(sig2).*eps;
            Pf_h0 = cumsum(u)+y0;
            Pf = Pf_h0;%cumsum(u)+y0;
    
            % Bubble and Burst component               
            % u1  = u(1:t1);        
            % B   = H \ ([B0; zeros(t1-1,1)]+  u1)-B0+y0; 
            % y   = det + [B; Pf(t1+1:end)];

            
            B0 = Pf_h0(t0);
            u1 = u(t0+1:t1);  
            if lam0==0
                B0=y0;
                u1 = u(1:t1);
            end
            
            
            
            
            B = HB \ ([rho*B0; zeros(n1-1,1)]+  u1);
            % Pf(t0+1:t1) =B;
            % y = det + Pf;


            w=rho.^(1:n1)';
            wi=1./w;  
            X=wi'*u1+B0;
                
            pos = X>0;%sign(X);
            if allow_neg_bubbles || t0==T
                pos = 1;
                %Pf(t0+1:t1) = B;
                %y = det + Pf;
            else      
                %pos = X>0;
                %Pf(t0+1:t1) =sign(X).*B-sign(B0).*2.*B0;
                %y = det + Pf;
                eps = randn(T,1);
            end

        end
        if lam0==0
            Pf(1:t1) = B;
        else
            Pf(t0+1:t1) = B;
        end

        %plot([Pf, Pf_h0])
        y = det + Pf;

   
end