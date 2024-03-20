############################ POISSON MASS PROBABILITY #########################



#Distributions can take many forms
#There is the NORMAL or GAUSSIAN distribution which is a BINOMIAL bell shaped curve
# Many types of data such as   sampling a population, flipping a coin,
# rolling dice where each sample is independent of its neighbor. This phenomenon
# is driven by what is called the "CENTRAL LIMIT THEROM". which says the more you test the closer to an
# an average we get...like taking 100 samples between -5 and 5 
# we can generate data of this type with stats::dnorm() from the STATS package

x= seq(-5,5, length = 100)
y = dnorm(x)
plot(x,y,xlab = "Sequence", ylab = "f(x) or dnorm(sequence)", main = "Normal Distribution")


# There is also the POISSON distribution
# Sometimes  data is based on time or maps to physical space. This type of data is often skewed to one side
# where most samples appear to the left of a mean. For example the state of a system over the course of 10 minutes
# or where an element exist on the genome. Data of this type is driven by the "POISSON PROBABILITY DENSITY FUNCTION"
# which says things of higher occurance probably occur earlier            P(x)= [(mu^x)* (e^x)]/x!
#    e = log constant 2.46
#    mu= an average of some kind (a base among all bases, read type among all reads)
#    x = a count of some kind (a base, a type of transcript, many states that can exist in time)
#    x! = x  factorial

# In sequencing this is used to assemble all the bases of a read into a cohesive sequence
# or to assemble transcripts. and reports back a probability on where these elements lie
# things being counted with high probability will exist earlier
# an average we get...like taking 100 samples the average of which is skewed around 10
# we can generate data of this type with stats::dnorm() from the STATS package
x<- 0:100
plot(dpois(x, 10), lwd = 2,
     main = "Poisson probability function",
     ylab = "P(X = x)", xlab = "Count")

plot(y,type = "h",xlab = "P(X=x)", ylab = "f(x) or dpois(sequence)", main = "Poisson Distribution")




# we can use the probablity mass function as a model and feed it a different type of average and different counts
# For example consider we have 4 types (counts) of read lengths.of a total of 10 reads. Theta is an avg of type and total 


type_1 (length= 800), (theta_1 = 4/10)
###################
###################
###################
###################

type_2 (length =100), (theta_2 3/10)
######
#####
#####


type_3 (length= 200), (theta_3 2/10)
#############
#############

type_4 (length= 1000), (theta_4 = 1/10)
###########################################



# the models would be where Y is the counts belonging to exon 1 2 3 and 4
Y_1 = poisson(L_1, theta_1)
Y_2 = poisson(L_2, theta_2)
Y_3 = poisson(L_3, theta_3)
Y_4 = poisson(L_4, theta_4)

# because exons are neighbors we can refine this a bit more with 
Y_1 = poisson(L_1 * theta_1 + L_1 * theta_2                                 )
Y_2 = poisson(L_2 * theta_1 + L_2 * theta_2 + L_2 * theta_3                 ) 
Y_3 = poisson(                L_3 * theta_2 + L_3 * theta_3 + L_2 * theta_4 )
Y_4 = poisson(                                L_4 * theta_3 + L_4 * theta_4 )
# the junctions would be similar
j1  = poisson(L_1 * theta_1 + L_1 * theta_2                                 )
j2  = poisson(L_2 * theta_1 + L_2 * theta_2 + L_2 * theta_3                 ) 
j3  = poisson(                L_3 * theta_2 + L_3 * theta_3 + L_2 * theta_4 )



# the above can be represented and modeled with the following matrix where Y/L = Poisson(M *theta) 

  {Y1/L1}          {1,1,0,0}
  {Y2/L2}          {1,1,1,0}  {theta_1}
  {Y3/L3}          {0,1,1,1}  {theta_2}
  {Y4/L4}= poisson {0,0,1,1} *{theta_3}
  {j1/L4}          {1,1,0,0}  {theta_4}
  {j2/L4}          {1,1,1,0}
  {j3/L4}          {0,1,1,1}


# HENCE using this model we can predict Y a vector of Poisson probabilities for exons and junctions if we have......
# L = a vector of lengths (L1,L2,L3,L4)
# theta = a vector of thetas (theta 1, theta 2,theta 3,theta 4)
# M = a matrix indicating which exons apply to a given theta 





l1 <- 100
l2 <- 200
l3 <- 300
l12 <- 100
l23 <- 100


lengths <- c(100, 200, 300, 100, 100)



mat <- cbind(c(1, 1, 0, 1, 0), 
             c(1, 1, 1, 1, 1), 
             c(0, 1, 1, 0, 1))


lengths %*% mat 



counts <- c(60, 320, 420, 60, 140)
w <- 1000



theta.hat <- c(1, 2, 3) / 10000


mat %*% theta.hat * lengths * w


LHS <- counts / (lengths * w)
lm.fit(mat, LHS)$coefficients
