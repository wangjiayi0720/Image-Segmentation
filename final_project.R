#####################################################################
# 2015 Spring JHU AMS 550.431 Final Project                         #
# Hao Wang                                                          #
# Jiayi Wang                                                        #
#####################################################################

setwd("/Users/haocabbage/Desktop/550.431")

# convert the training images to a 256x256 matrix
# which represents the observed x at each pixel
library(png)
training <- readPNG("training_ss_149.png")[, , 1] # training image
seg <- readPNG("training_seg_149.png")[, , 1] # segmented training image
test <- readPNG("test_ss_155.png")[, , 1] # test image
seg_test <- readPNG("test_seg_155.png")[, , 1] # segmented test image

rotate <- function(m) {
  # rotate the square matrix 90 degree to the right
  rev <- m
  for (i in 1:dim(m)[1]) {
    rev[i,] <- rev(m[,i])
  }
  return(rev)
}

training <- rotate(training)
seg <- rotate(seg)
test <- rotate(test)
seg_test <- rotate(seg_test)

# plot images 
par(mar=c(0,0,0,0))
image(training, col = gray(0:256/256), axes = FALSE)
image(seg, col = gray(0:256/256), axes = FALSE)
image(test, col = gray(0:256/256), axes = FALSE)
image(seg_test, col = gray(0:256/256), axes = FALSE)

# value transformation
trans <- function(m) {
  m_new <- as.vector(m)
  label_value <- sort(unique(m_new))
  for (i in 0:3) {
    m_new[m_new == label_value[i+1]] <- i
  }
  return(matrix(m_new, dim(m)[1], dim(m)[2]))  
}

seg <- trans(seg)
seg_test <- trans(seg_test)

# density of x in the original data
par(mar=c(3, 3, 2, 2))
hist(training, main = "Histogram of the training data")

#### ================================================================
#3.# 
####

# estimate mu_k, sigma_k via GDA
GDA.phi <- numeric(3)
GDA.mu <- numeric(3)
GDA.sigma <- numeric(3)

GDA.x <- training[31:180, 81:200]
GDA.y <- seg[31:180, 81:200]

# compute phi and mu
for (i in 1:dim(GDA.y)[1]) {
  for (j in 1:dim(GDA.y)[2]) {
    for (k in 1:3) {
      if (GDA.y[i,j] == k) {
        GDA.mu[k] <- GDA.mu[k] + GDA.x[i,j]
      }
    }
  }  
}

for (k in 1:3) {
  GDA.phi[k] <- length(GDA.y[GDA.y == k])/length(GDA.y[GDA.y != 0])
  GDA.mu[k] <- GDA.mu[k]/length(GDA.y[GDA.y == k])
}

# compute sigma
for (i in 1:dim(GDA.y)[1]) {
  for (j in 1:dim(GDA.y)[2]) {
    for (k in 1:3) {
      if (GDA.y[i,j] == k) {
        GDA.sigma[k] <- GDA.sigma[k] + (GDA.x[i,j] - GDA.mu[k])^2
      }
    }
  }  
}

for (k in 1:3) {
  GDA.sigma[k] <- sqrt(GDA.sigma[k]/length(GDA.y[GDA.y == k]))
}

# some helper functions
gibbs_sampling <- function(m, x, lambda, gamma, sweep) {
  # sample a graph with the given Gibbs parameters
  # m: pixel matrix
  # x: observed matrix
  # lambda, gamma: Gibbs parameters, 1x3 vectors
  # sweep: control the times of sweeping over the graph
  # sweep over the input matrix
  m_new <- m
  
  # expand gamma for easier reference
  gamma.m <- matrix(0, 3, 3)
  gamma.m[1,2] <- gamma[1]
  gamma.m[2,1] <- gamma[1]
  gamma.m[1,3] <- gamma[2]
  gamma.m[3,1] <- gamma[2]
  gamma.m[2,3] <- gamma[3]
  gamma.m[3,2] <- gamma[3]
  
  # sweeping
  for (sw in 1:sweep) {
    for (i in 1:dim(m_new)[1]) {
      for (j in 1:dim(m_new)[2]) {
        # start of y_s sampling ----------------------
        # check for black pixels
        if (m_new[i,j] != 0) {
          # get qualified neighbors
          neighbor <- c(m_new[i-1, j], m_new[i+1, j], 
                        m_new[i, j-1], m_new[i, j+1])
          neighbor <- neighbor[neighbor != 0]
          # get conditional prob. for y_s
          pr <- numeric(3)
          for (k in 1:3) {
            # compute psi function
            psi <- 0
            if (length(neighbor) > 0) {
              for (item in neighbor) {
                psi <- psi + gamma.m[k, item]
              }
            }
            pr[k] <- dnorm(x[i,j], GDA.mu[k], GDA.sigma[k])*
                     exp(lambda[k] + psi)
          }
          pr <- pr/sum(pr)
          #pr[is.na(pr)] <- 1
          # sample y_s
          m_new[i,j] <- sample(c(1,2,3), 1, prob=pr)
        }         
        # end of y_s sampling ------------------------        
      }
    }
  }
  return(m_new)
}

gamma_grad <- function(m) {
  # compute gradient part for updating gamma 
  # m: a pixel matrix
  toReturn <- c(0, 0, 0)
  for (i in 1:dim(m)[1]) {
    for (j in 1:dim(m)[2]) {
      if (m[i,j] != 0) {
        # get qualified neighbors
        neighbor <- c(m[i+1, j], m[i, j+1]) 
        neighbor <- neighbor[neighbor != 0]
        # compute gradient part
        if(length(neighbor) > 0) {
          for (item in neighbor) {
            temp <- sort(c(item, m[i,j]))
            if (identical(temp, c(1,2))) {
              toReturn[1] <- toReturn[1] + 1
            }
            if (identical(temp, c(1,3))) {
              toReturn[2] <- toReturn[2] + 1
            }
            if (identical(temp, c(2,3))) {
              toReturn[3] <- toReturn[3] + 1
            }
          }
        }
      }
    }
  }
  return(toReturn)
}

stoch_grad <- function(m, x, lambda, gamma, alpha, iter, sweep) {
  # estimate parameters using stochastic gradient
  # m: pixel matrix
  # x: observed matrix
  # lambda, gamma: Gibbs parameters, 1x3 vectors  
  # alpha: learning rate
  # iter: number of iterations
  # sweep: sweep times of Gibbs sampling
  # compute expected derivatives based on the true Gibbs dist.
  Ep0_lambda <- numeric(0)
  Ep0_theta <- numeric(0)
  for (k in 1:3) {
    Ep0_lambda[k] <- length(m[m==k])    
  }
  Ep0_theta <- gamma_grad(m)
  
  cat("Ep0_lambda:", Ep0_lambda, '\n', sep=" ")
  cat("Ep0_theta:", Ep0_theta, '\n', sep=" ")
  
  # construct a matrix to store parameters
  toReturn <- matrix(NA, 1, 6)
  toReturn[1, ] <- c(lambda, gamma)
  
  # construct a random matrix for gibbs sampling
  m_rand <- m
  for (i in 1:dim(m)[1]) {
    for (j in 1:dim(m)[2]) {
      if (m_rand[i,j] != 0) {
        m_rand[i,j] <- sample(c(1,2,3), 1, prob=c(1/3, 1/3, 1/3))
      }    
    }
  }

  # loop designated times 
  for (iteration in 1:iter) {
    # use the parameters from the previous turn
    lamb <- toReturn[iteration, ][1:3]
    gam <- toReturn[iteration, ][4:6]  
    
    # generate a matrix ~ Gibbs dist. 
    m_temp <- gibbs_sampling(m_rand, x, lamb, gam, sweep)
    
    # update parameters
    Ep_lambda <- numeric(0)
    Ep_theta <- numeric(0)
    for (k in 1:3) {
      Ep_lambda[k] <- length(m_temp[m_temp==k])    
    }
    Ep_theta <- gamma_grad(m_temp)
    
    lamb <- lamb + alpha*(Ep0_lambda - Ep_lambda)
    gam <- gam + alpha*(Ep0_theta - Ep_theta)
  
    # store parameters
    toReturn <- rbind(toReturn, c(lamb, gam))
    cat("iteration", iteration, ":", toReturn[iteration+1,], '\n', sep=" ")
    cat("Ep_lambda:", Ep_lambda, '\n', sep=" ")
    cat("Ep_theta:", Ep_theta, '\n', sep=" ")
  }
  return(toReturn)
}

error_rate <- function (m1, m2) {
  # compute the percentage of pixels that are mislabeled
  missed <- length(m2[m2 != m1])
  return(missed/length(m2[m2 != 0]))
}

# train without excessive black pixels
seg_crop <- seg[31:180, 81:200]
par(mar=c(0,0,0,0))
image(seg_crop, col = gray(0:256/256), axes = FALSE)

# select lambda and gamma by trial and error
m_trial <- gibbs_sampling(seg_crop, GDA.x, c(1,2,2), c(-3,-3,-4), 100)
par(mar=c(0,0,0,0))
image(m_trial, col = gray(0:256/256), axes = FALSE)
error_rate(seg_crop, m_trial)

# select lambda and gamma using stochastic gradient
# initialize Gibbs parameters
# gamma.0[1] = gamma_{1,2}
# gamma.0[2] = gamma_{1,3}
# gamma.0[3] = gamma_{2,3}
lambda.0 <- c(1, 2, 2)
gamma.0 <- c(-3, -3, -4) 

parameter.sg <- stoch_grad(seg_crop, GDA.x, lambda.0, gamma.0, 
                           0.001, 20, 40)

# simulate y with the optimal Gibbs parameters
lambda.sg <- c(1.4, 1.92, 1.71)
gamma.sg <- c(-0.6, -3, -0.6)
m_sg <- gibbs_sampling(seg_crop, GDA.x, lambda.sg, gamma.sg, 100)
par(mar=c(0,0,0,0))
image(m_sg, col = gray(0:256/256), axes = FALSE)
error_rate(seg_crop, m_sg)

#### ================================================================
#4.# 
####

# crop the training image for comparison
training_crop <- training[31:180, 81:200]
par(mar=c(0,0,0,0))
image(training_crop, col = gray(0:256/256), axes = FALSE)

# simulate 5 observed images
x.sim <- list(5)
for (x in 1:5) {
  m_temp <- GDA.x
  for (i in 1:dim(m_temp)[1]) {
    for (j in 1:dim(m_temp)[2]) {
      if (m_sg[i,j] != 0) {
        m_temp[i,j] <- rnorm(1, GDA.mu[m_sg[i,j]], 
                             GDA.sigma[m_sg[i,j]])
      }  
    }  
  }
  x.sim[[x]] <- m_temp
  title <- paste("task_4_x_", x, ".png", sep="")
  png(filename=title, width=150, height=120)
  par(mar=c(0,0,0,0))
  image(m_temp, col = gray(0:256/256), axes = FALSE)
  dev.off()
}

#### ================================================================
#5.# 
####

gibbs_sampling_freq <- function(m, x, lambda, gamma, sweep) {
  # sample a graph with the given Gibbs parameters and return
  # a matrix with each pixel taking the value with the highest
  # frequency
  
  # m: pixel matrix
  # x: observed matrix
  # lambda, gamma: Gibbs parameters, 1x3 vectors
  # sweep: control the times of sweeping over the graph
  # sweep over the input matrix
  m_new <- m
  
  # construct a list of 3 matrices to store frequency for each 
  # value the pixel can take
  m_freq <- list()
  for (k in 1:3) {
    m_freq[[k]] <- matrix(0, dim(m)[1], dim(m)[2])  
  }
  
  # construct a 3x3 matrix of gamma for easier reference
  gamma.m <- matrix(0, 3, 3)
  gamma.m[1,2] <- gamma[1]
  gamma.m[2,1] <- gamma[1]
  gamma.m[1,3] <- gamma[2]
  gamma.m[3,1] <- gamma[2]
  gamma.m[2,3] <- gamma[3]
  gamma.m[3,2] <- gamma[3]
  
  # sweeping
  for (sw in 1:sweep) {
    for (i in 1:dim(m_new)[1]) {
      for (j in 1:dim(m_new)[2]) {
        # start of y_s sampling ----------------------
        # check for black pixels
        if (m_new[i,j] != 0) {
          # get qualified neighbors
          neighbor <- c(m_new[i-1, j], m_new[i+1, j], 
                        m_new[i, j-1], m_new[i, j+1])
          neighbor <- neighbor[neighbor != 0]
          # get conditional prob. for y_s
          pr <- numeric(3)
          for (k in 1:3) {
            # compute psi function
            psi <- 0
            if (length(neighbor) > 0) {
              for (item in neighbor) {
                psi <- psi + gamma.m[k, item]
              }
            }
            pr[k] <- dnorm(x[i,j], GDA.mu[k], GDA.sigma[k])*
              exp(lambda[k] + psi)
          }
          pr <- pr/sum(pr)
          #pr[is.na(pr)] <- 1
          # sample y_s
          m_new[i,j] <- sample(c(1,2,3), 1, prob=pr)
          # update frequency matrices
          freq <- m_new[i,j]
          m_freq[[freq]][i,j] <- m_freq[[freq]][i,j] + 1  
        }         
        # end of y_s sampling ------------------------        
      }
    }
  }
  # pick the value with the highest frequency for each pixel
  for (i in 1:dim(m_new)[1]) {
    for (j in 1:dim(m_new)[2]) {
      if (m_new[i,j] != 0) {
        hi_freq <- 0
        for (k in 1:3) {
          if(m_freq[[k]][i,j] > hi_freq) {
            m_new[i,j] <- k 
          }  
        } 
      }
    }
  }
  return(m_new)
}

# estimate y with simulated x
y.est <- list()
for (y in 1:5) {
  y.est[[y]] <- gibbs_sampling_freq(seg_crop, x.sim[[y]], 
                                    lambda.sg, gamma.sg, 100)  
  cat("y_", y, ": ", error_rate(seg_crop, y.est[[y]]), "\n", sep="")
  title <- paste("task_5_y_", y, ".png", sep="")
  png(filename=title, width=150, height=120)
  par(mar=c(0,0,0,0))
  image(y.est[[y]], col = gray(0:256/256), axes = FALSE)
  dev.off()
}

# estimate with gibbs_sampling instead of gibbs_sampling_freq
y.sg <- list()
for (y in 1:5) {
  y.sg[[y]] <- gibbs_sampling(seg_crop, x.sim[[y]], 
                                    lambda.sg, gamma.sg, 200)  
  cat("y_", y, ": ", error_rate(seg_crop, y.sg[[y]]), "\n", sep="")
  title <- paste("task_5_y_sg_", y, ".png", sep="")
  png(filename=title, width=150, height=120)
  par(mar=c(0,0,0,0))
  image(y.sg[[y]], col = gray(0:256/256), axes = FALSE)
  dev.off()
}

#### ================================================================
#6.# 
####

# Belief Propagation 
compute_msg <- function(i, j, neb, m, x, msg_box, lambda, gamma) {
  # compute the raw message from pixel (i,j) to the designated neighbor
  # i, j: the position of y_s
  # neb: define which neighbor {Down, Right, Up, Left} = {1,2,3,4}
  # m: the pixel matrix
  # x: the observed matrix
  # msg_box: the messege box
  # lambda, gamma: gibbs parameters
  # expand gamma for easier reference
  gamma.m <- matrix(0, 3, 3)
  gamma.m[1,2] <- gamma[1]
  gamma.m[2,1] <- gamma[1]
  gamma.m[1,3] <- gamma[2]
  gamma.m[3,1] <- gamma[2]
  gamma.m[2,3] <- gamma[3]
  gamma.m[3,2] <- gamma[3]
  
  msg <- numeric(3)  
  for (t in 1:3) {
    msg_part <- 0
    for (s in 1:3) { 
      # compute phi(x_s, y_s)  
      phi.ys <- exp(lambda[s]*dnorm(x[i,j], GDA.mu[s], GDA.sigma[s]))
      # compute psi(y_s, y_t)
      psi.ys <- exp(gamma.m[s,t])
      # gather messages from other neighbors
      neb_msg <- 1
      if (neb == 1) { # exclude down
        if (m[i,j+1] != 0) { # right
          neb_msg <- neb_msg*msg_box[[i]][[j+1]][[4]][s]
        }   
        if (m[i-1,j] != 0) { # up
          neb_msg <- neb_msg*msg_box[[i-1]][[j]][[1]][s] 
        }   
        if (m[i,j-1] != 0) { # left
          neb_msg <- neb_msg*msg_box[[i]][[j-1]][[2]][s] 
        }     
      } 
      if (neb == 2) { # exclude right
        if (m[i+1,j] != 0) { # down
          neb_msg <- neb_msg*msg_box[[i+1]][[j]][[3]][s]
        }   
        if (m[i-1,j] != 0) { # up
          neb_msg <- neb_msg*msg_box[[i-1]][[j]][[1]][s] 
        }   
        if (m[i,j-1] != 0) { # left
          neb_msg <- neb_msg*msg_box[[i]][[j-1]][[2]][s] 
        }     
      }
      if (neb == 3) { # exclude up
        if (m[i+1,j] != 0) { # down
          neb_msg <- neb_msg*msg_box[[i+1]][[j]][[3]][s]
        } 
        if (m[i,j+1] != 0) { # right
          neb_msg <- neb_msg*msg_box[[i]][[j+1]][[4]][s]
        }   
        if (m[i,j-1] != 0) { # left
          neb_msg <- neb_msg*msg_box[[i]][[j-1]][[2]][s] 
        }     
      }
      if (neb == 4) { # exclude left
        if (m[i+1,j] != 0) { # down
          neb_msg <- neb_msg*msg_box[[i+1]][[j]][[3]][s]
        }
        if (m[i,j+1] != 0) { # right
          neb_msg <- neb_msg*msg_box[[i]][[j+1]][[4]][s]
        }   
        if (m[i-1,j] != 0) { # up
          neb_msg <- neb_msg*msg_box[[i-1]][[j]][[1]][s] 
        }        
      } 
      # finish up
      msg_part <- msg_part + phi.ys*psi.ys*neb_msg  
    }
    msg[t] <- msg_part
  }
  #normalize message
  msg <- msg/sum(msg)
  return(msg)
}

BP_alg <- function(m, x, lambda, gamma, iter) {
  # implementation of Belief Propagation
  # m: a pixel matrix
  # x: observed matrix
  # lambda, gamma: gibbs parameters
  # iter: number of iterations
  
  # scan the graph and initialize all messages  
  msg_box <- list() # it's going to be a 3D structure
  # in which msg_box[[i]][[j]][[v]][k]
  # is msg_{(i,j) ->  neighbor v}(y_v = k)
  
  for (i in 1:dim(m)[1]) {
    # construct message box
    msg_box[[i]] <- list()
    
    for (j in 1:dim(m)[2]) {
      # continue constructing message box
      msg_box[[i]][[j]] <- list()
      for (v in 1:4) {
        msg_box[[i]][[j]][[v]] <- c()
      }
      # check for validity of pixels
      if (m[i,j] != 0) {
        # detect neighbors and initialize msg 
        # direction: down -> right -> up -> left
        if (m[i+1,j] != 0) {
          msg_box[[i]][[j]][[1]] <- c(1/3, 1/3, 1/3)
        } 
        if (m[i,j+1] != 0) {
          msg_box[[i]][[j]][[2]] <- c(1/3, 1/3, 1/3)
        }   
        if (m[i-1,j] != 0) {
          msg_box[[i]][[j]][[3]] <- c(1/3, 1/3, 1/3)
        }   
        if (m[i,j-1] != 0) {
          msg_box[[i]][[j]][[4]] <- c(1/3, 1/3, 1/3)
        }   
      }
    }
  }
  
  # iterate designated times
  for (iteration in 1:iter) {
    # update all messages
    for (i in 1:dim(m)[1]) {
      for (j in 1:dim(m)[2]) {
        # check for validity of pixels
        if (m[i,j] != 0) {
          # pass new messages to neighbors
          # direction: down -> right -> up -> left
          if (m[i+1,j] != 0) {
            msg_box[[i]][[j]][[1]] <- 
              compute_msg(i, j, 1, m, x, msg_box, lambda, gamma) 
          } 
          if (m[i,j+1] != 0) {
            msg_box[[i]][[j]][[2]] <- 
              compute_msg(i, j, 2, m, x, msg_box, lambda, gamma)   
          }   
          if (m[i-1,j] != 0) {
            msg_box[[i]][[j]][[3]] <- 
              compute_msg(i, j, 3, m, x, msg_box, lambda, gamma)   
          }   
          if (m[i,j-1] != 0) {
            msg_box[[i]][[j]][[4]] <- 
              compute_msg(i, j, 4, m, x, msg_box, lambda, gamma) 
          }   
        }
      }
    }  
  }
  # return the message information
  return(msg_box)
}

compute_belief <- function(i, j, m, x, msg_box, lambda) {
  # compute the belief vector for y_s
  # i, j: the position of y_s
  # m: the pixel matrix
  # x: the observed matrix
  # msg_box: the messege box
  # lambda: gibbs parameter
  
  belief <- numeric(3)
  for (k in 1:3) {
    # compute phi(x_s, y_s)
    phi.ys <- exp(lambda[k]*dnorm(x[i,j], GDA.mu[k], GDA.sigma[k]))
    # gather messages from neighbors
    neb_msg <- 1
    if (m[i+1,j] != 0) { # down
      neb_msg <- neb_msg*msg_box[[i+1]][[j]][[3]][k]
    }  
    if (m[i,j+1] != 0) { # right
      neb_msg <- neb_msg*msg_box[[i]][[j+1]][[4]][k]
    }   
    if (m[i-1,j] != 0) { # up
      neb_msg <- neb_msg*msg_box[[i-1]][[j]][[1]][k] 
    }   
    if (m[i,j-1] != 0) { # left
      neb_msg <- neb_msg*msg_box[[i]][[j-1]][[2]][k] 
    }
    belief[k] <- phi.ys*neb_msg  
  }
  return(belief)
}

BP_seg <- function(m, x, msg_box, lambda) {
  # segmentation based on BP algorithm
  m_new <- m
  for (i in 1:dim(m_new)[1]) {
    for (j in 1:dim(m_new)[2]) {
      # check vadility of pixels
      if (m_new[i,j] != 0) {
        # get belief of y_s
        belief <- compute_belief(i, j, m_new, x, msg_box, lambda)
        # maximize p(y_s=k|x)
        m_new[i,j] <- which.max(belief)
      }
    }
  }
  return(m_new)
}

# apply BP algorithm using the training data
training.msg <- BP_alg(seg_crop, GDA.x, lambda.sg, gamma.sg, 20)

# segment with BP message
m_BP <- BP_seg(seg_crop, GDA.x, training.msg, lambda.sg)
par(mar=c(0,0,0,0))
image(m_BP, col = gray(0:256/256), axes = FALSE)
error_rate(m_BP, seg_crop)

# more comparisons 
# visualize which pixels are different
comp_vis <- seg_crop
comp_vis[comp_vis == m_BP] <- 0
par(mar=c(0,0,0,0))
image(comp_vis, col = gray(0:256/256), axes = FALSE)

# 1(y_s_BP = 1 & y_s = 1)
count1 <- 0
# 1(y_s_BP = 2 & y_s = 2)
count2 <- 0
# 1(y_s_BP = 3 & y_s = 3)
count3 <- 0

for (i in 1:dim(seg_crop)[1]) {
  for (j in 1:dim(seg_crop)[2]){
    if(seg_crop[i,j] == 1 && m_BP[i,j] == 1) {
      count1 <- count1 + 1
    }
    if(seg_crop[i,j] == 2 && m_BP[i,j] == 2) {
      count2 <- count2 + 1
    }
    if(seg_crop[i,j] == 3 && m_BP[i,j] == 3) {
      count3 <- count3 + 1
    }
  }
}

# P(y_s_BP = 1|y_s = 1)
cat("P(y_s_BP = 1|y_s = 1) = ", count1/length(seg_crop[seg_crop == 1]), "\n")
# P(y_s_BP = 2|y_s = 2)
cat("P(y_s_BP = 2|y_s = 2) = ", count2/length(seg_crop[seg_crop == 2]), "\n")
# P(y_s_BP = 3|y_s = 3)
cat("P(y_s_BP = 3|y_s = 3) = ", count3/length(seg_crop[seg_crop == 3]), "\n")

#### ================================================================
#7.# 
####

# crop the test images for segmentation 
test_crop <- test[31:180, 81:200]
par(mar=c(0,0,0,0))
image(test_crop, col = gray(0:256/256), axes = FALSE)

seg_test_crop <- seg_test[31:180, 81:200]
par(mar=c(0,0,0,0))
image(seg_test_crop, col = gray(0:256/256), axes = FALSE)

# apply gibbs sampler to the test image
m_test <- gibbs_sampling(seg_test_crop, test_crop, lambda.0, gamma.0, 200)
par(mar=c(0,0,0,0))
image(m_test, col = gray(0:256/256), axes = FALSE)
cat("y_test:", error_rate(seg_test_crop, m_test), "\n")

# apply BP algorithm to the test image
test.msg <- BP_alg(seg_test_crop, test_crop, lambda.sg, gamma.sg, 20)
m_BP_test <- BP_seg(seg_test_crop, test_crop, test.msg, lambda.sg)
par(mar=c(0,0,0,0))
image(m_BP_test, col = gray(0:256/256), axes = FALSE)
error_rate(m_BP_test, seg_test_crop)

# more comparisons 
# visualize which pixels are different
comp_vis_test <- seg_test_crop
comp_vis_test[comp_vis_test == m_BP_test] <- 0
par(mar=c(0,0,0,0))
image(comp_vis_test, col = gray(0:256/256), axes = FALSE)

# 1(y_s_BP_test = 1 & y_s = 1)
count1_test <- 0
# 1(y_s_BP_test = 2 & y_s = 2)
count2_test <- 0
# 1(y_s_BP_test = 3 & y_s = 3)
count3_test <- 0

for (i in 1:dim(seg_test_crop)[1]) {
  for (j in 1:dim(seg_test_crop)[2]){
    if(seg_test_crop[i,j] == 1 && m_BP_test[i,j] == 1) {
      count1_test <- count1_test + 1
    }
    if(seg_test_crop[i,j] == 2 && m_BP_test[i,j] == 2) {
      count2_test <- count2_test + 1
    }
    if(seg_test_crop[i,j] == 3 && m_BP_test[i,j] == 3) {
      count3_test <- count3_test + 1
    }
  }
}

# P(y_s_BP_test = 1|y_s = 1)
cat("P(y_s_BP_test = 1|y_s = 1) = ", count1_test/length(seg_test_crop[seg_test_crop == 1]), "\n")
# P(y_s_BP_test = 2|y_s = 2)
cat("P(y_s_BP_test = 2|y_s = 2) = ", count2_test/length(seg_test_crop[seg_test_crop == 2]), "\n")
# P(y_s_BP_test = 3|y_s = 3)
cat("P(y_s_BP_test = 3|y_s = 3) = ", count3_test/length(seg_test_crop[seg_test_crop == 3]), "\n")

################### =================================================
# EM Segmentation # 
###################

EM_alg <- function(v, m, p, s, tol) {
  # implementation of EM algorithm
  # v: data vector
  # m: initial means
  # p: initial proportion of means
  # s: initial standard deviation
  # tol: tolerance level for convergence
  n <- length(v) # number of data points
  k <- length(m) # number of labels
  Q_value <- 0 # expected log likelihood
  
  # loop until convergence
  while(TRUE) {
    # E-step
    gamma <- matrix(NA, n, k)
    for (j in 1:k) {
      gamma[,j] <- 1/(s[j]*sqrt(2*pi))*exp(-(v-m[j])^2/(2*s[j]^2))*p[j]
    }
    gamma <- gamma/rowSums(gamma)
    
    # update Q_value
    temp <- matrix(NA, n, k)
    for (j in 1:k) {
      temp[,j] <- (log(1/(s[j]*sqrt(2*pi))*exp(-(v-m[j])^2/(2*s[j]^2)))
                   +log(p[j]))*gamma[,j]
    }
    Q_value_new <- sum(temp) 
    cat(paste("Q Value = ", Q_value, "\n", sep=""))
    if (abs(Q_value - Q_value_new) <= tol) { break }
    else { Q_value <- Q_value_new }
    
    # M-step
    m_old <- m # use m in the last iteration
    p <- colSums(gamma)/n
    s <- sqrt((colSums(gamma*(replicate(k, v)-t(replicate(n, m_old)))^2))/
                (colSums(gamma)))
    m <- (colSums(gamma*replicate(k, v)))/(colSums(gamma))
  }
  
  return(list(mean=m, prop=p, sigma=s))
}

# train without the black pixels
black <- training[1]
dat <- as.vector(training)
dat_nb <- dat[dat != black] # nb stands for no black

# apply EM alg to the training data
result <- EM_alg(dat_nb, sort(runif(3)), c(1/3, 1/3, 1/3), 
                 c(0.2, 0.2, 0.2), 0.0001)

#EM parameters
mu <- result$mean 
phi <- result$prop
sigma <- result$sigma

# generate y with the x and EM parameters
# assign 0 to black pixels and {1, 2, 3} to other
# pixels according to the posterior dist. of y_s
m_em <- dat
for (i in 1:length(m_em)) {
  if (m_em[i] == black) {
    m_em[i] <- 0
  } else { 
    # compute posterior dist.
    pr <- numeric(3)
    for (j in 1:length(mu)) {
      pr[j] <- dnorm(m_em[i], mu[j], sigma[j])*phi[j]  
    }
    pr <- pr/sum(pr)
    # sample y_s
    m_em[i] <- sample(c(1,2,3), 1, prob=pr)
  }
}

m_em <- matrix(m_em, 256, 256)
par(mar=c(0,0,0,0))
image(m_em, col = gray(0:256/256), axes = FALSE)

# crop it for visual comparison
m_em_crop <- m_em[31:180, 81:200]
par(mar=c(0,0,0,0))
image(m_em_crop, col = gray(0:256/256), axes = FALSE)

# trim the matrix to find error rate
for (i in 1:256) {
  for (j in 1:256) {
    if (m_em[i,j] != 0 && seg[i,j] == 0) {
      m_em[i,j] <- 0
    }
  }
}
error_rate(m_em_crop, seg_crop)
