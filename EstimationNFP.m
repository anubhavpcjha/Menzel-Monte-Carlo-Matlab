clear;

options = optimset('Algorithm','sqp');

rand('seed',1234);
%rand('seed',5678);

B = 20;
N=[100,200,500,1000,2000];

thNFP=ones(8,B,length(N));
thMPEC=ones(8,B,length(N));

timeNFP=ones(B,length(N));
timeMPEC=ones(B,length(N));

i=1;
for n = N
    b=1;
    while b<B+1
        % first regressor is constant, second is group indicator
        n2 = n/2;
        X = repmat([1 0;1 1],n2,1);
        Z = repmat([1 0;1 1],n2,1);

        
        
        beta_w = [0.5;0.5;0;0;0;0;0;1];
        beta_m = [0.5;0.5;0;0;0;0;0;1];
       
        % new formulation of outside option
        U_star = makeIndex(X,Z,beta_w);
        V_star = makeIndex(Z,X,beta_m)';
        
        
        eta = -log(-log(rand(n,n)));
        zeta = -log(-log(rand(n,n)));
        
        J = round(sqrt(n));
        eta0 = repmat(max(-log(-log(rand(J,n))))',1,n);
        zeta0 = repmat(max(-log(-log(rand(J,n)))),n,1);
        
        eta = eta - eta0;
        zeta = zeta - zeta0;

       
        
        U = U_star + eta;
        V = V_star + zeta;
        

        mu = Gale_Shapley(U,V);

        % Gamma is 4 by 1 vector, first two components are women (0 and 1),
        % last two components are men (0 and 1)
        
        k = length(beta_w);
        
        theta_0 = [beta_w;beta_m;zeros(4,1)];
      
        loglik = @(theta)loglik3(theta,mu,X,Z);
        loglik(theta_0)
        
        n
        b
        
        % NFP
         A_res = eye(2*k+4); 
        A_res(k+1,k+1) = 0;
        A_res(1,k+1)= - 1;
        A_res(2*k,2*k) = 0;
        A_res(k,2*k)= - 1; 
        A_res(k+2,k+2) = 0;
        A_res(2,k+2)= - 1;
        A_res(2*k+1:2*k+4,2*k+1:2*k+4) = zeros(4,4);
        
        b_res = zeros(2*k+4,1);
        
        
        r=cputime;        
        options = optimset('Algorithm','active-set','MaxIter',200);
        [th_hat,~,~] = fmincon(loglik,theta_0,[],[],A_res,b_res,[],[],[],options);
        timeNFP(b,i)=cputime-r;
        
        thNFP(:,b,i)=th_hat(1:8,1);
        
        
        %MPEC
        A_ineq = zeros(2*k+4);
        A_ineq(2*k+1:2*k+4,2*k+1:2*k+4) = -eye(4);
        

       
             
        
        
        loglik0 = @(theta)loglik2(theta,mu,X,Z);
        psi_diff = @(theta)psi_fp1(theta,X,Z);
        
        t=cputime;
        options = optimset('Algorithm','active-set','MaxIter',200);
        [th_hat,~,~] = fmincon(loglik0,theta_0,A_ineq,b_res,A_res,b_res,[],[],psi_diff,options);
        timeMPEC(b,i)=cputime-t;
        
        if norm(th_hat(1:8,1))<50
            thMPEC(:,b,i)=th_hat(1:8,1);
            b=b+1;
        end
       
     
        
        
        
    end
    
    
    
    
    i=i+1;
    
end
    nfp=mean(thNFP,2);
    mpec=mean(thMPEC,2);
    
    meantimeNFP=mean(timeNFP,1);
    meantimeMPEC=mean(timeMPEC,1);
   
    nfpresid=thNFP-repmat(beta_w,1,B,length(N));
    mpecresid=thMPEC-repmat(beta_w,1,B,length(N));
    
    stderrNFP=mean(nfpresid.^2,2);
    stderrMPEC=mean(mpecresid.^2,2);
    

        