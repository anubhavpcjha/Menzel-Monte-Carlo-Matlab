clear;

options = optimset('Algorithm','sqp');

rand('seed',1234);
%rand('seed',5678);

i=1;

n = 200;






Alph_mat = 0.5;

    

    B = 10;
    
    beta_sim = zeros(B,9);
    
    tic
   

        
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

%         eta = randn(n,n);
%         zeta = randn(n,n);
%         eta0 = repmat(max(randn(J,n))',1,n);
%         zeta0 = repmat(max(randn(J,n)),n,1);
%  %       sig_n = sqrt(lambertw(n/2/pi));
%         b_n = norminv(1 - 1/sqrt(n));
%         sig_n = sqrt(n)*normpdf(b_n); 
%         eta = sig_n*(eta - eta0);
%         zeta = sig_n*(zeta - zeta0);
        
        
        U = U_star + eta;
        V = V_star + zeta;
        

        mu = Gale_Shapley(U,V);

        % Gamma is 4 by 1 vector, first two components are women (0 and 1),
        % last two components are men (0 and 1)
        
        k = length(beta_w);
        
        theta_0 = [beta_w;beta_m;zeros(4,1)];
       %theta_0 = [zeros(2*k,1);ones(4,1)];
        loglik = @(theta)loglik3(theta,mu,X,Z);
        loglik(theta_0)
        
        A_ineq = zeros(2*k+4);
        A_ineq(2*k+1:2*k+4,2*k+1:2*k+4) = -eye(4);
        

        A_res = eye(2*k+4); 
        A_res(k+1,k+1) = 0;
        A_res(1,k+1)= - 1;
        A_res(2*k,2*k) = 0;
        A_res(k,2*k)= - 1; 
        A_res(k+2,k+2) = 0;
        A_res(2,k+2)= - 1;
        A_res(2*k+1:2*k+4,2*k+1:2*k+4) = zeros(4,4);

        b_res = zeros(2*k+4,1);
        
        tic
        options = optimset('Algorithm','active-set','MaxIter',200);
        [th_hat,fval,exitflag] = fmincon(loglik,theta_0,[],[],A_res,b_res,[],[],[],options);
        toc
        
        
        loglik0 = @(theta)loglik2(theta,mu,X,Z);
        psi_diff = @(theta)psi_fp1(theta,X,Z);
        
        tic
        options = optimset('Algorithm','active-set','MaxIter',200);
        [th_hat,fval,exitflag] = fmincon(loglik0,theta_0,A_ineq,b_res,A_res,b_res,[],[],psi_diff,options);
        toc
        
        