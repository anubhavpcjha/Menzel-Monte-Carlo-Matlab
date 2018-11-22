function mu = Gale_Shapley(U,V)

% first argument proposing side
% proposing side: individuals correspond to rows

[nw,nm] = size(U);
U_temp = U;
nmax = 10*nw*nm;

for i = 1:nmax
    
 % check dimensions!   
 % no need to worry about ties
    Prop = (U_temp==repmat(max(max(U_temp,[],2),0),1,nm));
    Rej = (Prop.*V < Prop.*max(0,repmat(max(Prop.*V,[],1),nw,1)));
    U_temp(Rej&Prop) = -1;
 
    if sum(sum(Rej))==0
        break
    end;

end;

mu = (U_temp==repmat(max(max(U_temp,[],2),0),1,nm));

end