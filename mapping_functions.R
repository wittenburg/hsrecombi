### general mapping functions under the assumption of no Chromatid Interference (i.e. 0<= theta <= 0.5) 

# Haldane (1919)
haldane=function(x, inverse = F){
  y=c()
  if(inverse){ # Morgan -> theta
    for (i in 1:length(x)){
      y[i]=0.5*(1-exp(-2*x[i]))
    }  
  } else{ # theta -> Morgan
    for (i in 1:length(x)){
      theta=min(x[i], 0.499)
      y[i]=-0.5*log(1-2*theta)
    }
  }
  y
}

# Kosambi (1944)
kosambi=function(x, inverse = F){
  y=c()
  if(inverse){ # Morgan -> theta
    for (i in 1:length(x)){
      y[i]=0.5*tanh(2*x[i])
    }
  } else{ # theta -> Morgan
    for (i in 1:length(x)){
      theta=min(x[i], 0.499)
      y[i]=1/4*log((1+2*theta)/(1-2*theta))
    }
  }
  y
}


# Rao (1977); # p=0 would match to Morgan map function (p=0.25 to Carter, p=0.5 to Kosambi, p=1 to Haldane)
rao=function(p,x, inverse = F){
  y=c()
  if(inverse){ # Morgan -> theta
    y <- rep(NA, length(x))
  } else{ # theta -> Morgan
    for(i in 1:length(x)){
      theta=min(x[i], 0.499)
      y[i]=(p*(2*p-1)*(1-4*p)*log(1-2*theta)+16*p*(p-1)*(2*p-1)*atan(2*theta)+
              2*p*(1-p)*(8*p+2)*atanh(2*theta)+6*(1-p)*(1-2*p)*(1-4*p)*theta)/6
    }
  }
  y
}


# approximation to inverse of Rao's mapping function
rao.inv <- function(p, x){ # Morgan -> theta
  theta <- c()
  for(i in 1:length(x)){
    opt <- optim(par = 0.1, fn = function(th){(x[i] - rao(p, th))^2}, method = 'Brent', lower = 0, upper = 0.49, control = list(reltol = 1e-5)) 
    theta[i] <- opt$par
  }
  theta
}


# Felsenstein (1979)
felsenstein=function(K,x, inverse = F){
  y=c()
  if(inverse){ # Morgan -> theta
    for(i in 1:length(x)){
      y[i]=(1-exp(2*(K-2)*x[i]))/2/(1-(K-1)*exp(2*(K-2)*x[i]))
    }  
  } else{ # theta -> Morgan
    for(i in 1:length(x)){
      theta=min(x[i],0.499)
      y[i]=1/2/(K-2)*log((1-2*theta)/(1-2*(K-1)*theta))
    }
  }
  y
}


# Karlin's Binomial (1984); N=1 Morgan's mapping function
karlin <- function(N,x, inverse = F){
  y <- c()
  if(inverse){ # Morgan -> theta
    for (i in 1:length(x)){
      y[i] <- ifelse(x[i] < N / 2, 0.5 * (1 - (1 - 2 * x[i] / N)^N), 1 / 2)
    }  
  } else{ # theta -> Morgan
    for (i in 1:length(x)){
      theta <- min(x[i], 0.499)
      y[i] <- 0.5 * N * (1 - (1 - 2 * theta)^(1 / N))
    }
  }
  y  
}


# Carter & Falconer (1951)
carter <- function(x, inverse = F){
  y <- c()
  if(inverse){ # Morgan -> theta
    y <- rep(NA, length(x))
  } else{ # theta -> Morgan
    for (i in 1:length(x)){
      theta <- min(x[i], 0.499)
      y[i] <- 0.25 * (atan(2 * theta) + atanh(2* theta))
    }
  }
  y  
}
