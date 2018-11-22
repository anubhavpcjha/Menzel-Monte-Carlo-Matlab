function LL = loglik3(theta,mu,X,Z)

    nw = length(X);
    nm = length(Z);
    k = length(theta)/2-2;
    beta_w = theta(1:k,1);
    beta_m = theta(k+1:2*k,1);
    Gamma_old = theta(2*k+1:2*k+4,1);
    
    
    iter=1;
    maxiter=1000;
    tol=10^(-8);
    while norm(psi_fp4(theta,Gamma_old,X,Z)-Gamma_old)>tol && iter<=maxiter
        Gamma_old=psi_fp4(theta,Gamma_old,X,Z);
        iter=iter+1;
    end
    Gamma_old=psi_fp4(theta,Gamma_old,X,Z);
    
    

    Gam_w_mat = log(ones(nw,1) + Gamma_old(1,1)*(1-X(:,2)) + Gamma_old(2,1)*X(:,2));
    Gam_m_mat = log(ones(nm,1) + Gamma_old(3,1)*(1-Z(:,2)) + Gamma_old(4,1)*Z(:,2));
    
    
%     U_star = makeIndex(X,Z,beta_w);
%     V_star = makeIndex(X,Z,beta_m);

    U_star = makeIndex(X,Z,beta_w);
    V_star = makeIndex(Z,X,beta_m)';
   
    W_star = mu.*(U_star + V_star);
    
    mar_w = sum(mu,2);
    mar_m = sum(mu,1)';
    
    LL = -(2*sum(sum(W_star)) - sum((1+mar_w).*Gam_w_mat) - sum((1+mar_m).*Gam_m_mat));
    
end
