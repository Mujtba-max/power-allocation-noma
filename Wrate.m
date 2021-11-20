% Please update this for the complex v !!!

function R = Wrate(NC, NU, H, v, alpha, nvar)
R = zeros(NU, NC);
for c = 1:NC
  for i = 1:NU
      inter = 0;
      for k=1:NC
          if k~=c
              temp = H(i,c,k)^2*sum(v(:, k).^2);
              inter = inter + temp;
          end
      end
      if i == 1
          R(i, c) = alpha(i, c) * log2(1+(H(i, c, c)*v(i, c))^2/(nvar+inter));
      else
          R(i, c) = alpha(i, c) * log2(1+(H(i, c, c)*v(i, c))^2/(inter+sum((H(1:i-1 , c, c).*v(1:i-1, c)).^2)+nvar));
      end
  end
end

end