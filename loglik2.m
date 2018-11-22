function LL = loglik2(theta,mu,X,Z)

    nw = length(X);
    nm = length(Z);
    k = length(theta)/2-2;
    beta_w = theta(1:k,1);
    beta_m = theta(k+1:2*k,1);
    Gamma = theta(2*k+1:2*k+4,1);
    
    

    Gam_w_mat = log(ones(nw,1) + Gamma(1,1)*(1-X(:,2)) + Gamma(2,1)*X(:,2));
    Gam_m_mat = log(ones(nm,1) + Gamma(3,1)*(1-Z(:,2)) + Gamma(4,1)*Z(:,2));
    
    
%     U_star = makeIndex(X,Z,beta_w);
%     V_star = makeIndex(X,Z,beta_m);

    U_star = makeIndex(X,Z,beta_w);
    V_star = makeIndex(Z,X,beta_m)';
   
    W_star = mu.*(U_star + V_star);
    
    mar_w = sum(mu,2);
    mar_m = sum(mu,1)';
    
    LL = -(2*sum(sum(W_star)) - sum((1+mar_w).*Gam_w_mat) - sum((1+mar_m).*Gam_m_mat));
    
end

