function R = alphaTdmaRates(NC,NU,h,alpha,Pmax,nvar)
R = zeros(NU,NC);
g = abs(h).^2;
nv = nvar/NU;
Pt = Pmax/NU;
for c=1:NC
    for u=1:NU
        SINR = (g(u,c,c)*Pt)/(sum(g(u,c,:).*Pt)+nv);
        R(u,c) = alpha(u,c)*log2(1+SINR);
    end
end