function Pn_theta = n_coverage_probability_monte_carlo_SIC(theta, n, SIC)
  ## Pn_theta: Joint pdf of the SIR and successive interference cancellation
  ## Input:
  ## theta: SIR threshold of a successful transmission
  ## n: transmitter's index
  ## SIC: set '0' to evaluate the n-probability, otherwise evaluate SIC-SIR

  K = max(n, 2); # Interference cancellation stages (K >= n)  
  M = 5e03; # M: number of Monte Carlo samples
  tic # Measure and print the evaluation time
  nukappa = 2; # "Tildekappa x nu" -parameter
  tau = min(theta, 0.5); # Signal detection threshold (tau <= theta)
  tau_prime = tau / (1 + tau);
  
  ## Integration lower bound
  theta_prime = theta / (1 + theta);
  
  ## Upper bound for the series expansion of the joint PDF
  imax = max(0,ceil(1 / theta_prime - n - 1));  
  
  ## If n==1 and SIC == 0 evaluate the built-in numerical integration, 
  ## which is much more efficient compared to the Monte Carlo method
  ## in the simple coverage region theta >=1
  if n == 1 && theta >= 1 && SIC == 0 && true
    mu_prime = @(t) nukappa ./ t .* (1-t) .^ (nukappa - 1);    
    Pn_theta = integral(mu_prime, theta_prime, 1);

    ## Else perform Monte Carlo integration over z'
  else
    Pn_theta = 0;
    if SIC == 0;
      Pn_theta = deltasic(n);
    else
      for k = n : K
	Pn_theta = Pn_theta + deltasic(k);
      endfor
    endif    
  endif
  toc

  function res = deltasic(k)
    ## Initialize an array to store the calculated function values
    deltasic_k_values = zeros(1, M);
    for iii = 1 : M
      if mod(iii, 250) == 0
	iii
      end
      
      ## The conditioning and the integration region are defined separately
      ## without and with SIC 
      if SIC == 0
	## Generate n random samples for z' from the uniform distribution on
	## (theta_prime, 1) and define the integration domain
	z_prime_samples = theta_prime + (1 - theta_prime) * rand(1, k);
	volume_of_domain = (1 - theta_prime) ^ k; #Integration domain volume
	is_sorted = all(diff(z_prime_samples) < 0); # Sorting condition	
	cond = is_sorted;
      else
	
	## Generate n random samples for z' from the uniform distribution
	z_prime_samples = rand(1, k);
	volume_of_domain = 1; #Integration domain volume
	is_sorted = all(diff(z_prime_samples) < 0); # Sorting condition
	## Condition 1
	cumulative_sum = cumsum(z_prime_samples);
	shifted_cumulative_sum = [0, cumulative_sum(1:end-1)];
	condition_vector = z_prime_samples...
			   + tau_prime .* shifted_cumulative_sum ...
			   > tau_prime;
	cond1 = all(condition_vector);
	## Condition 2
	sum_part = sum(z_prime_samples(1:k-1)) - z_prime_samples(n);
	cond2 = (k>n) * (z_prime_samples(n)...
		       + theta_prime * sum_part < theta_prime) + (k==n);
	## Condition 3
	sum_part = sum(z_prime_samples(1:k)) - z_prime_samples(n);
	cond3 = z_prime_samples(n) + theta_prime * sum_part > theta_prime;
	cond = is_sorted * cond1 * cond2 * cond3;
      endif

      ## Calculate the joint PDF f'_(n) for the current sample             
      ## as the series expansion of the joint PDF
      deltasic_value = 0;      
      for j = 0 : imax
	if sum(z_prime_samples) + j * z_prime_samples(end) < 1 && cond
	  mu_prime_val = calculate_mu_prime(z_prime_samples, j);
	else
	  mu_prime_val = 0;
	endif
	term = ((-1) ^ j / factorial(j)) * mu_prime_val;
	deltasic_value = deltasic_value + term;
      end      
      deltasic_k_values(iii) = deltasic_value;
    end

    ## Apply the Monte Carlo integration formula
    average_deltasic_k = mean(deltasic_k_values);
    res = volume_of_domain * average_deltasic_k;
    
  endfunction
  ## A nested function for the mu_prime
  function mu_prime_val = calculate_mu_prime(z_prime_vector, j)
    z_prime_dim=length(z_prime_vector);
    ## z_prime_dim: The dimension of the z' vector provided.
    ## i: The number of variables (zeta) to integrate over.
    ## nukappa: The parameter of the gamma process.
    ## z_prime_vector: The vector of n z' variables (z'_1, ..., z'_n).
    ## M: The number of Monte Carlo samples.
    

    ## The integration domain for each zeta_j is [z'_n, 1].
    z_prime_n = z_prime_vector(end);
    
    ## Calculate the volume of the integration domain
    volume = (1 - z_prime_n) ^ j;
    ## Initialize a vector to store the values of the integrand
    integrand_values = zeros(1, M);
    
    ## Perform M Monte Carlo samples
    for index = 1 : M
      
      ## Generate ordered uniform r.v.'s in rectangle
      zeta_samples = z_prime_n + (1 - z_prime_n) * rand(1, j);
      zeta_samples=sort(zeta_samples,'descend');
      full_t_vector = [z_prime_vector, zeta_samples];
      
      ## Check if within the simplex.
      sum_is_le_1 = sum(full_t_vector) <= 1;
      
      if sum_is_le_1
	## If conditions are met, calculate the value of the density
	product_term = prod(full_t_vector .^ (-1));
	sum_term = sum(full_t_vector);
	density_val = (nukappa) ^ (z_prime_dim + j) * product_term ...
		      * (1 - sum_term) ^ (nukappa - 1);
	integrand_values(index) = density_val;
      else
	## If conditions are not met, the density is 0
	integrand_values(index) = 0;
      end
    end
    
    ## Calculate the final integral value using the Monte Carlo formula
    average_integrand = mean(integrand_values);
    mu_prime_val = average_integrand * volume;
  endfunction
endfunction



