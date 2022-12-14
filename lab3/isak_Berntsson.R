# By Jose M. Pe?a and Joel Oskarsson.
# For teaching purposes.
# jose.m.pena@liu.se.

#####################################################################################################
# Q-learning
#####################################################################################################

#install.packages("ggplot2")
#install.packages("vctrs")
library(ggplot2)

# --------------------- section A -------------------
arrows <- c("^", ">", "v", "<")
action_deltas <- list(c(1,0), # up
                      c(0,1), # right
                      c(-1,0), # down
                      c(0,-1)) # left

vis_environment <- function(iterations=0, epsilon = 0.5, alpha = 0.1, gamma = 0.95, beta = 0){
  
  # Visualize an environment with rewards. 
  # Q-values for all actions are displayed on the edges of each tile.
  # The (greedy) policy for each state is also displayed.
  # 
  # Args:
  #   iterations, epsilon, alpha, gamma, beta (optional): for the figure title.
  #   reward_map (global variable): a HxW array containing the reward given at each state.
  #   q_table (global variable): a HxWx4 array containing Q-values for each state-action pair.
  #   H, W (global variables): environment dimensions.
  
  df <- expand.grid(x=1:H,y=1:W)
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,q_table[x,y,1],NA),df$x,df$y)
  df$val1 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,q_table[x,y,2],NA),df$x,df$y)
  df$val2 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,q_table[x,y,3],NA),df$x,df$y)
  df$val3 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,q_table[x,y,4],NA),df$x,df$y)
  df$val4 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) 
    ifelse(reward_map[x,y] == 0,arrows[GreedyPolicy(x,y)],reward_map[x,y]),df$x,df$y)
  df$val5 <- as.vector(foo)
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,max(q_table[x,y,]),
                                     ifelse(reward_map[x,y]<0,NA,reward_map[x,y])),df$x,df$y)
  df$val6 <- as.vector(foo)
  
  print(ggplot(df,aes(x = y,y = x)) +
          scale_fill_gradient(low = "white", high = "green", na.value = "red", name = "") +
          geom_tile(aes(fill=val6)) +
          geom_text(aes(label = val1),size = 4,nudge_y = .35,na.rm = TRUE) +
          geom_text(aes(label = val2),size = 4,nudge_x = .35,na.rm = TRUE) +
          geom_text(aes(label = val3),size = 4,nudge_y = -.35,na.rm = TRUE) +
          geom_text(aes(label = val4),size = 4,nudge_x = -.35,na.rm = TRUE) +
          geom_text(aes(label = val5),size = 10) +
          geom_tile(fill = 'transparent', colour = 'black') + 
          ggtitle(paste("Q-table after ",iterations," iterations\n",
                        "(epsilon = ",epsilon,", alpha = ",alpha,"gamma = ",gamma,", beta = ",beta,")")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_x_continuous(breaks = c(1:W),labels = c(1:W)) +
          scale_y_continuous(breaks = c(1:H),labels = c(1:H)))
  
}

GreedyPolicy <- function(x, y){
  
  # Get a greedy action for state (x,y) from q_table.
  #
  # Args:
  #   x, y: state coordinates.
  #   q_table (global variable): a HxWx4 array containing Q-values for each state-action pair.
  # 
  # Returns:
  #   An action, i.e. integer in {1,2,3,4}.
  
  return(sample(which(q_table[x,y,]==max(q_table[x,y,])),1))
  
}

EpsilonGreedyPolicy <- function(x, y, epsilon){
  
  # Get an epsilon-greedy action for state (x,y) from q_table.
  #
  # Args:
  #   x, y: state coordinates.
  #   epsilon: probability of acting randomly.
  # 
  # Returns:
  #   An action, i.e. integer in {1,2,3,4}.
  
  # Your code here.
  if (runif(1) > epsilon){ # P(random) = 1-epsilon
    return(sample(1:4,1))
  }else{

    return(sample(which(q_table[x,y,]==max(q_table[x,y,])),1))
    
  }
  
  
}

transition_model <- function(x, y, action, beta){
  
  # Computes the new state after given action is taken. The agent will follow the action 
  # with probability (1-beta) and slip to the right or left with probability beta/2 each.
  # 
  # Args:
  #   x, y: state coordinates.
  #   action: which action the agent takes (in {1,2,3,4}).
  #   beta: probability of the agent slipping to the side when trying to move.
  #   H, W (global variables): environment dimensions.
  # 
  # Returns:
  #   The new state after the action has been taken.
  
  delta <- sample(-1:1, size = 1, prob = c(0.5*beta,1-beta,0.5*beta))
  final_action <- ((action + delta + 3) %% 4) + 1
  foo <- c(x,y) + unlist(action_deltas[final_action])
  foo <- pmax(c(1,1),pmin(foo,c(H,W)))
  
  #print(foo)
  return (foo)
}

q_learning <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95, 
                       beta = 0){
  
  # Perform one episode of Q-learning. The agent should move around in the 
  # environment using the given transition model and update the Q-table.
  # The episode ends when the agent reaches a terminal state.
  # 
  # Args:
  #   start_state: array with two entries, describing the starting position of the agent.
  #   epsilon (optional): probability of acting greedily.
  #   alpha (optional): learning rate.
  #   gamma (optional): discount factor.
  #   beta (optional): slipping factor.
  #   reward_map (global variable): a HxW array containing the reward given at each state.
  #   q_table (global variable): a HxWx4 array containing Q-values for each state-action pair.
  # 
  # Returns:
  #   reward: reward received in the episode.
  #   correction: sum of the temporal difference correction terms over the episode.
  #   q_table (global variable): Recall that R passes arguments by value. So, q_table being
  #   a global variable can be modified with the superassigment operator <<-.
  
  # Your code here.
  
  
  state = start_state
  
  correction = 0
  
  episode_correction = 0
  repeat{

    action = EpsilonGreedyPolicy(state[1],state[2],epsilon = epsilon )
    #print(action)
    
    new_state = transition_model(state[1],state[2],action = action, beta = beta)
    
    #print(new_state)
    
    reward = reward_map[new_state[1], new_state[2]]
    #print(reward)
    
    Q_SA = q_table[state[1],state[2],action]
    
    correction = alpha*(reward + gamma*max(q_table[new_state[1], new_state[2],]) - Q_SA)
    
    episode_correction = episode_correction + correction/alpha
    
    
    q_table[state[1],state[2],action] <<- Q_SA + correction     
    #q_table[state[1],state[2],action] <<- Q_SA + alpha*(reward + gamma*max(q_table[new_state[1], new_state[2],]) - Q_SA)
        
    if(reward!=0){
      # End episode.
      return (c(reward,episode_correction))
    }
    
    state = new_state  
  }
  
}

#####################################################################################################
# Q-Learning Environments
#####################################################################################################

# Environment A (learning)
set.seed(1337)
H <- 5
W <- 7

reward_map <- matrix(0, nrow = H, ncol = W)
reward_map[3,6] <- 10
reward_map[2:4,3] <- -1


q_table <- array(0,dim = c(H,W,4))


vis_environment()

for(i in 1:10000){
  foo <- q_learning(start_state = c(3,1))
  
  if(any(i==c(10,100,1000,10000,100000)))
    vis_environment(i)
}

"After 10 episodes the agent has learned that some of the squares close to the 'wall' are likely to lead to a negative result."

"No, state[1,3] has greedy policy to go straight into the 'wall'."

"Yes somewhat, the lower path around the wall has lower values. We could implement epsiloan annealing,
thus forcing the agent to expolre he space in the earlier iterations. Before finding the optimal paths "

""
#---------------------------------------- envB Learning --------------------------

set.seed(1337)
# Environment B (the effect of epsilon and gamma)

H <- 7
W <- 8

reward_map <- matrix(0, nrow = H, ncol = W)
reward_map[1,] <- -1
reward_map[7,] <- -1
reward_map[4,5] <- 5
reward_map[4,8] <- 10

q_table <- array(0,dim = c(H,W,4))

vis_environment()

MovingAverage <- function(x, n){
  
  cx <- c(0,cumsum(x))
  rsum <- (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]) / n
  
  return (rsum)
}
#set.seed(1337)

for(j in c(0.5,0.75,0.95)){
  q_table <- array(0,dim = c(H,W,4))
  reward <- NULL
  correction <- NULL
  
  for(i in 1:30000){
    foo <- q_learning(epsilon=5, gamma = j, start_state = c(4,1))
    reward <- c(reward,foo[1])
    correction <- c(correction,foo[2])
  }
  
  vis_environment(i, gamma = j)
  title = paste(" Epsilon = ", .5, "\n",
                " gamma = ", j , "\n",
                " Reward"
                )
  plot(MovingAverage(reward,100),type = "l", main=title)
  
  title = paste(" Epsilon = ", .5, "\n",
                " gamma = ", j , "\n",
                " Correction"
  )
  plot(MovingAverage(correction,100),type = "l",main=title)
}

for(j in c(0.5,0.75,0.95)){
  q_table <- array(0,dim = c(H,W,4))
  reward <- NULL
  correction <- NULL
  
  for(i in 1:30000){
    foo <- q_learning(epsilon = 0.1, gamma = j, start_state = c(4,1))
    reward <- c(reward,foo[1])
    correction <- c(correction,foo[2])
  }
  
  vis_environment(i, epsilon = 0.1, gamma = j)
  
  title = paste(" Epsilon = ", .1, "\n",
                " gamma = ", j, "\n",
                "Reward"
  )
  plot(MovingAverage(reward,100),type = "l",main=title)
  
  title = paste(" Epsilon = ", .1, "\n",
                " gamma = ", j, "\n",
                "Correction"
  )
  plot(MovingAverage(correction,100),type = "l",main=title)
  
}

"A low gamma penalizes the agent based on the length of the path. Since 10 is further away from the start than 5. A lower gamma means the agent will focus on finding the path to 5. 
We also see that the value of the squares is much lower further away from the rewards when gamma is low"

"A higher epsilon makes the algorithm converge faster, e = .5 -> ~3000it.  e = .1 -> ~15000it "


#------------------------------- env c --------------------------------


# Environment C (the effect of beta).

H <- 3
W <- 6

reward_map <- matrix(0, nrow = H, ncol = W)
reward_map[1,2:5] <- -1
reward_map[1,6] <- 10

q_table <- array(0,dim = c(H,W,4))

vis_environment()

for(j in c(0,0.2,0.4,0.66)){-
  q_table <- array(0,dim = c(H,W,4))
  
  for(i in 1:10000)
    foo <- q_learning(gamma = 0.6, beta = j, start_state = c(1,1))
  
  vis_environment(i, gamma = 0.6, beta = j)
}

"A high beta can be interpreted as the result of our actions being uncertain. This means we are incentivised to take a longer route, far away for the negative rewards when beta is high. 
In this run 1-gamma is supstantial and there is some tension between wanting to take a short path and wanting to take a safe path. In general, a high beta will take us on a longer path,m
further from the negative rewrds"






#----------------------------------REINFORCE-------------------------------------

#----------------------------------- Environemtn D---------------------------

" The agent does learn policies that lead to the goal. However, it is obvious that it doesn't understand the problem.
The probabilities of going into the goal from an adjacent square are sometimes very clsoe to zero and the agent may heavily favor
taking one route over another even though they are equal in length"

"For the Q-learning algorithm to handle new different environments the state has to be described in relative terms I.E. we are (dX, dY) squares from the goal, 
rather than: we are in square (x,y). This solution would, I think, work very well. Although that problem is simple enough you could easily just calculate an optimal solution"


#-------------------------- Environment E --------------------------------
"The agent shows a very clear bias for moving upwards. This can obviously be attributed to the mismatch in training and validation data. 
The agent is simply not tested on the same type of data it has been trained on. "

" Again, when we train on such small datasets, and they show such clear differences the results are boudn to be very different. 
The first model shows no clear bias in what action to take,m wheras he second one does.
"

