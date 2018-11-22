function Xb = makeIndex(Xwomen,Xmen,theta_w)

[nw,k] = size(Xwomen);
[nm,k] = size(Xmen);

tw_own = theta_w(1:k,1);
tw_oth = theta_w(k+1:2*k,1);
tw_diff = theta_w(2*k+1:3*k,1);
tw_dist = theta_w(3*k+1:4*k,1);

% need to check dimensions/reshape
Xb = reshape((repmat(Xwomen,nm,1)-reshape(repmat(Xmen',nw,1),k,nm*nw)').^2*tw_dist,nw,nm);
%Xb = reshape((repmat(Xwomen,nm,1)-kron(Xmen,ones(nw,1))).^2*tw_dist,nw,nm);



Xb = Xb + repmat(Xwomen*tw_own,1,nm) + repmat(tw_oth',nw,1)*Xmen';
Xb = Xb + repmat(Xwomen*tw_diff,1,nm) - repmat(tw_diff',nw,1)*Xmen';

end