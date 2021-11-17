function res = suff(user1,user2,Pmax,h,NC,cell,nvar)
h = abs(h);
hcell = h(cell,cell,:);
nvar1 = nvar(cell,user1);
nvar2 = nvar(cell,user2);
Q = zeros(1,NC);
for i=1:NC
    if(i~=cell)
        Q(1,i) = (hcell(user1)/hcell(user2)) < (h(i,cell,user1)/h(i,cell,user2)) + 0;
    end
end
RHS = sum(Q(1,:).*(Pmax(1,:)./(nvar1*nvar2)).*(h(:,cell,user1)'*hcell(user2)-h(:,cell,user2)'*hcell(user1)));
LHS = hcell(user1)/nvar1-hcell(user2)/nvar2;
res = (LHS>RHS) + 0;
end