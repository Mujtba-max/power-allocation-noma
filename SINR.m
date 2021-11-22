function S = SINR(NC, NU, H, v, nvar)
  % I think the dimentions here are old and reversed ! Does this function still in use?
S = zeros(NC, NU);
for c = 1:NC
  for i = 1:NU
      inter = 0;
      for k=1:NC
          if k~=c
              temp = H(i,c,k)^2*sum(v(k,:).^2);
              inter = inter + temp;
          end
      end
      if i == 1
          S(c, i) = (H(i, c, c)*v(c, i))^2/(nvar+inter);
      else
          S(c, i) = (H(i, c, c)*v(c, i))^2/(inter+sum((H(1:i-1 , c, c).*v(c, 1:i-1)').^2)+nvar);
      end
  end
end
end