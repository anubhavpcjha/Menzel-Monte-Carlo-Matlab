function psi_eq = psi_fp4(theta,Gamma_old,X,Z)


    nw = length(X);
    nm = length(Z);
    k = length(theta)/2-2;
    beta_w = theta(1:k,1);
    beta_m = theta(k+1:2*k,1);
    Gamma = Gamma_old;

    Gam_w_mat = ones(nw,1) + Gamma(1,1)*(1-X(:,2)) + Gamma(2,1)*X(:,2);
    Gam_m_mat = ones(nm,1) + Gamma(3,1)*(1-Z(:,2)) + Gamma(4,1)*Z(:,2);
    
    
    % need to set x values to zero and one, respectively for FP equation,
    % only need nm times number of types.
    
    Xw = [1 0;1 1];
    
%     U_star = makeIndex(Xw,Z,beta_w);
%     V_star = makeIndex(Xw,Z,beta_m);
%     W_star_w = U_star + V_star;
    
    U_star = makeIndex(Xw,Z,beta_w);
    V_star = makeIndex(Z,Xw,beta_m)';
    W_star_w = U_star + V_star;
    
    Zm = [1 0;1 1];
    
%     U_star = makeIndex(X,Zm,beta_w);
%     V_star = makeIndex(X,Zm,beta_m);
%     W_star_m = U_star + V_star;

    U_star = makeIndex(X,Zm,beta_w);
    V_star = makeIndex(Zm,X,beta_m)';
    W_star_m = U_star + V_star;

    
    Psi_w = mean(exp(W_star_w')./repmat(Gam_m_mat,1,2));
    Psi_m = mean(exp(W_star_m)./repmat(Gam_w_mat,1,2));
    
   
    psi_eq = [Psi_w,Psi_m]';
        
end

