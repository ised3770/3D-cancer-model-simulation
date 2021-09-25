function [Time_exp, Lya_exp] = lyapunov(n, func, ode_45, time_start, step_size, time_end, l_start, out);

% Lyapunov calculation
% A. Wolf, J. B. Swift, H. L. Swinney, and J. A. Vastano,
% "Determining Lyapunov Exponents from a Time Series".

n1 = n; 
n2 = n1*(n1+1);

% Number steps
step = round((time_end-time_start)/step_size);
 
% Allocation
y = zeros(n2,1); 
mem = zeros(n1,1); 
y0 = y;
gsc = mem;
znorm = mem;

y(1:n) = l_start(:);
for i = 1:n1 y((n1+1)*i) = 1.0; end;
% Start time
t = time_start;

% Main loop
for lya_index = 1:step 
  [T,Y] = feval(ode_45,func,[t t+step_size],y);  
  
  t = t+step_size;
  y = Y(size(Y,1),:);
  for i = 1:n1 
      for j = 1:n1 y0(n1*i+j)=y(n1*j+i); end;
  end;

% Orthonormal basis, Gram-Schmidt:
  znorm(1) = 0.0;
  for j=1:n1 znorm(1) = znorm(1)+y0(n1*j+1)^2; end;
  znorm(1) = sqrt(znorm(1));
  for j = 1:n1 y0(n1*j+1) = y0(n1*j+1)/znorm(1); end;
  for j = 2:n1
      for k = 1:(j-1)
          gsc(k) = 0.0;
          for l = 1:n1 gsc(k) = gsc(k) + y0(n1*l+j) * y0(n1*l+k); end;
      end;
 
      for k = 1:n1
          for l = 1:(j-1)
              y0(n1*k+j) = y0(n1*k+j) - gsc(l)*y0(n1*k+l);
          end;
      end;
      znorm(j) = 0.0;
      for k = 1:n1 znorm(j) = znorm(j) + y0(n1*k+j)^2; end;
      znorm(j) = sqrt(znorm(j));
      for k=1:n1 y0(n1*k+j) = y0(n1*k+j) / znorm(j); end;
  end;
  
% Update
  for k=1:n1 mem(k) = mem(k) + log(znorm(k)); end;
  
%Normalization of exponent
  for k=1:n1 
      lp(k) = mem(k)/(t-time_start); 
  end;
  
% Output
  if lya_index == 1
     Lya_exp = lp;
     Time_exp = t;
  else
     Lya_exp = [Lya_exp; lp];
     Time_exp = [Time_exp; t];
  end;
  
  % Print result
  if (mod(lya_index,out == 0))
     fprintf('t=%6.4f',t);
     for k=1:n1 fprintf(' %10.6f',lp(k)); end;
     fprintf('\n');
  end;
  for i = 1:n1 
      for j = 1:n1
          y(n1*j+i) = y0(n1*i+j);
      end;
  end;
end;