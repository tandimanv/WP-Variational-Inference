#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(rstanarm)
library(ggplot2)
library(extraDistr)
library(tidyverse)
library(gridExtra)
library(DT)
library(viridis)
library(tictoc)


# Define UI for application that draws a histogram
ui <- fluidPage(
      
      # Application title
      titlePanel("Variational Inference (VI) for Unit-Variance Univariate GMM with K=3"),
      # slider inputs at the top for initial m and s2 values
      fluidRow(
            column(3,wellPanel(
                  helpText("\\(x_i|c_i,\\boldsymbol{\\mu} \\sim N(c_i^T \\boldsymbol{\\mu},1), \\, i = 1,...,n\\)"),
                  helpText("\\(\\mu_k \\sim N(0,\\sigma^2), \\, k=1,...,K\\)"),
                  helpText("\\(c_i \\sim \\mbox{Categorical} (\\frac{1}{K},...,\\frac{1}{K}), \\, i = 1,...,n\\)"),
                  helpText("\\( \\,\\)"),
                  helpText("\\(\\mbox{Mean-Field VI}\\)"),
                  helpText("\\(q(\\boldsymbol{\\mu,c}) = \\prod_{k=1}^K {q(\\mu_k;m_k, s_k^2)} \\prod_{i=1}^n {q(c_i;\\psi_i)}\\)"),
                  helpText("\\(\\mbox{using CAVI updates}\\)"),
                  br(),
                  fluidRow(
                        splitLayout(
                              numericInput("n", label = "Sample size:",
                                           value = 100, min = 30, max = 500,
                                           step = 10),
                              numericInput("threshold", label = "Threshold:",
                                           value = 0.005, min = -Inf, max = Inf,
                                           step = 0.01)
                        )
                  ),
                  fluidRow(
                        sliderInput("sig2", label = withMathJax("\\(\\sigma^2\\):"),
                                    min = 1, max = 20, value = 10, step = 1)
                  )
            )
            ),
            column(3, wellPanel(
                  sliderInput("m1", label = withMathJax("Initial \\(m_1\\):"),
                              min = -20, max = 20, value = -5, step = 1),
                  sliderInput("m2", label = "Initial \\(m_2\\):",
                              min = -20, max = 20, value = 0, step = 1),
                  sliderInput("m3", label = "Initial \\(m_3\\):",
                              min = -20, max = 20, value = 10, step = 1)
            )
            ),
            column(3, wellPanel(
                  sliderInput("s2_1", label = "Initial \\(s_1^2\\):",
                              min = 0, max = 20, value = 3, step = 1),
                  sliderInput("s2_2", label = "Initial \\(s_2^2\\):",
                              min = 0, max = 20, value = 8, step = 1),
                  sliderInput("s2_3", label = "Initial \\(s_3^2\\):",
                              min = 0, max = 20, value = 2, step = 1)
            )
            ),
            column(3, wellPanel(
                  sliderInput("psi1", label = "Initial \\(\\psi_1\\):",
                              min = 0, max = 1, value = 0.2, step = 0.05),
                  uiOutput("psi2"),
                  uiOutput("psi3"),
                  br()
            ) 
            )
      ),
      # True data generating values
      fluidRow(
            column(3, wellPanel(
                  h4("True Data Generating Values"),
                  sliderInput("mu1", label = "\\(\\mu_1\\):",
                              min = -10, max = 10, value = -5, step = 1),
                  sliderInput("mu2", label = "\\(\\mu_2\\):",
                              min = -10, max = 10, value = 0, step = 1),
                  sliderInput("mu3", label = "\\(\\mu_3\\):",
                              min = -10, max = 10, value = 3, step = 1)
            ) 
            ),
            # The output
            column(8,
                   tabsetPanel(type = "tabs",
                               tabPanel("Plot", plotOutput("distPlot")),
                               tabPanel("Progress", plotOutput("process")),
                               tabPanel("ELBO", plotOutput("elbo")),
                               tabPanel("Table", 
                                        column(8,DT::dataTableOutput("table")))
                   ))
      )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
      
      k <- 3
      output$psi2 <- renderUI({
            sliderInput("psi2", label = withMathJax("Initial \\(\\psi_2\\):"),
                        min = 0, max = 1 - input$psi1, value = 0.3, step = 0.05)
      })
      
      output$psi3 <- renderUI({
            # paste("Initial \\(\\psi_3\\) =", 1 - (input$psi1 + input$psi2), 
            #       sep = " ")
            withMathJax(paste("Initial \\(\\psi_3\\) :", 
                              1 - (input$psi1 + input$psi2),
                              sep = " "))
      })
      
      
      # Create reactive expression
      # generating the data and predictions
      output_dfsim <- reactive({
            # fix sigma^2 (hyperparameter for mu_k), k, and n
            sig2 <- input$sig2
            k <- 3
            n <- input$n
            # generate data based on the input mus and ps
            mu <- c(input$mu1, input$mu2, input$mu3)
            mu <- mu[order(mu)]
            # set cluster assignment
            set.seed(30)
            clust <- rcat(n,c(1/3, 1/3, 1/3))
            df_sim <- tibble(data = rep(0,n),
                             cluster = factor(clust))
            # set data per cluster assignment
            for (i in 1:n){
                  set.seed(29+i)
                  df_sim[i,1] <- rnorm(1,mu[clust[i]],1)
            }
            truth <- list(df_sim = df_sim,
                          mu = mu,
                          sig2 = sig2)
      })
      
      output_combo <- reactive({
            truth <- output_dfsim()
            df_sim <- truth$df_sim
            sig2 <- truth$sig2
            k <- 3
            n <- input$n
            
            
            # update and elbo functions
            update_psi <- function(m,s2,x){
                  psi <- matrix(0,nrow = length(x), ncol = length(m))
                  for (i in (1:length(x))){
                        psi[i,] <- exp(x[i]*m-(s2+m^2)/2)
                        psi[i,] <- psi[i,]/sum(psi[i,])
                  }
                  return(psi)
            } 
            
            update_s2 <- function(sig2,psi){
                  s2 <- 1/(1/sig2 + colSums(psi))
                  return(s2)
            }
            
            update_m <- function(sig2,psi,x){
                  s2 <- 1/(1/sig2 + colSums(psi))
                  m <- colSums(psi*x)*s2
                  return(m)
            }
            
            elbo <- function(m,s2,psi, x){
                  n<- length(x)
                  k <- length(m)
                  a <- sum(-.5*log(2*pi*sig2)-1/(2*sig2)*(s2+m^2))
                  for_b <- rep(0,n)
                  for (i in (1:n)){
                        for_b[i] <- -0.5*log(2*pi)-0.5*x[i]^2+x[i]*sum(psi[i,]*m)-.5*sum(psi[i,]*(s2+m^2))
                  }
                  b <- -n*log(k) + sum(for_b)
                  c <- sum(apply(psi*log(psi),1, sum))
                  d <- sum(-0.5*log(2*pi*s2)-0.5)
                  return(a+b-c-d)
            }
            
            # storage matrices
            threshold <- input$threshold
            max_iter <- 500
            psi_iter <- array(0, c(n,k,max_iter))
            # psi_iter <- matrix(0,nrow = n, ncol = k*max_iter) #every iter is a nx3 input
            m_iter <- matrix(0,nrow = max_iter, ncol = k)
            s2_iter <- matrix(0,nrow = max_iter, ncol = k)
            elbo_iter <- rep(0, max_iter)
            # note that if the loop doesn't reach max_iter, you'll
            # have some 0s at the remaining rows or elements
            
            m_inits <- c(input$m1,input$m2, input$m3)
            s2_inits <- c(input$s2_1,input$s2_2, input$s2_3)
            psi_inits <- matrix(rep(c(input$psi1, input$psi2, 1-(input$p1+input$p2)),n), 
                                nrow = n, ncol = k, byrow = T)
            
            # check ELBO for inits
            elbo_iter[1] <- elbo(m=m_inits,
                                 s2=s2_inits, 
                                 psi = psi_inits, 
                                 x= df_sim$data) 
            
            
            m_iter[1,] <- m_inits
            s2_iter[1,] <- s2_inits
            psi_iter[,,1] <- psi_inits
            
            
            i <- 2
            
            tic("Runtime")
            psi_iter[,,i] <- update_psi(m = m_iter[i-1,],
                                        s2 = s2_iter[i-1,],
                                        x = df_sim$data)
            s2_iter[i,] <- update_s2(sig2 = sig2, 
                                     psi =  psi_iter[,,i])
            m_iter[i,] <- update_m(sig2 = sig2, psi = psi_iter[,,i],
                                   x = df_sim$data)
            elbo_iter[i] <- elbo(m= m_iter[i,],
                                 s2 = s2_iter[i,],
                                 psi = psi_iter[,,i],
                                 x = df_sim$data)
            
            while((elbo_iter[i]-elbo_iter[i-1]) > threshold & i < max_iter){
                  i = i+1
                  psi_iter[,,i] <- update_psi(m = m_iter[i-1,],
                                              s2 = s2_iter[i-1,],
                                              x = df_sim$data)
                  s2_iter[i,] <- update_s2(sig2 = sig2, 
                                           psi =  psi_iter[,,i])
                  m_iter[i,] <- update_m(sig2 = sig2, psi = psi_iter[,,i],
                                         x = df_sim$data)
                  elbo_iter[i] <- elbo(m= m_iter[i,],
                                       s2 = s2_iter[i,],
                                       psi = psi_iter[,,i],
                                       x = df_sim$data)
            }
            toc(log = TRUE, quiet = TRUE)
            log.txt <- tic.log(format = TRUE)
            tic.clearlog()
            
            # re-order the groups from small to large
            urut <- order(m_iter[i,])
            m_iter <- m_iter[,urut]
            s2_iter <- s2_iter[,urut]
            for(s in (1:i)){
                  hold <- psi_iter[,,s]
                  psi_iter[,,s] <- hold[,urut]
            }
            combo <- list(
                  elbo_iter = elbo_iter,
                  psi_iter = psi_iter,
                  s2_iter = s2_iter,
                  m_iter = m_iter,
                  i = i,
                  runtime = unlist(log.txt)
            )
            combo
      }) 
      
      output$distPlot <- renderPlot({
            combo <- output_combo()
            elbo_iter <- combo$elbo_iter
            psi_iter <- combo$psi_iter
            s2_iter <- combo$s2_iter
            m_iter <- combo$m_iter
            psi <- combo$psi_iter
            i <- combo$i
            runtime <- combo$runtime
            
            truth <- output_dfsim()
            df_sim <- truth$df_sim
            mu <- truth$mu
            df_sim <- df_sim %>% 
                  mutate(data = round(data, 5),
                         predicted = as.factor(apply(psi[,,i],1,which.max)))
            # visualize Final
            
            # order rows of df_sim by data values
            idx <- order(df_sim$data)
            df_sim2 <- df_sim[idx,]
            
            binwidth <- diff(range(df_sim2$data))/35  
            
            ggplot(df_sim2, aes(x = data, fill = predicted)) +
                  geom_dotplot(dotsize = 0.5, binwidth = binwidth,
                               aes(color =  cluster, stroke = 3, alpha = 0.5)) +
                  theme_bw() + scale_color_viridis(discrete = TRUE, option = "E")+
                  scale_fill_viridis(discrete = TRUE) + 
                  scale_y_continuous(NULL, breaks = NULL) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[i,1], sd = sqrt(s2_iter[i,1])),
                                colour = 'grey', size = 1, alpha = 0.7, n = 20000) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[i,2], sd = sqrt(s2_iter[i,2])),
                                colour = 'grey', size = 1, alpha = 0.7, n = 20000) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[i,3], sd = sqrt(s2_iter[i,3])),
                                colour = 'grey', size = 1, alpha = 0.7, n = 20000) +
                  geom_vline(xintercept = mu, color = viridis(3, option = "E"), 
                             # linetype = 'dotted',
                             size = 1) + guides(alpha = FALSE) +
                  labs(fill = "Predicted \nAssignment", 
                       color = "True \nAssignment") + xlab("") + ylab("") +
                  # labs(caption = expression(paste("Total Iterations:", i, "        ",
                  #                                                      "Classification Error Rate:",
                  #                                                      100*mean(df_sim$cluster !=  as.numeric(df_sim$predicted)),
                  #                                                      "%"))) +
                  # theme(plot.caption = element_text(face = "bold", hjust = 0.5,
                  #                                   size = rel(1.5)))
                  labs(caption = bquote(atop(bold("Total Iterations: "~.(i)~"         Assignment Error Rate:"~
                                                        .(100*mean(df_sim$cluster !=  as.numeric(df_sim$predicted)))~
                                                        "%     "~.(runtime)),"Optimized Approximations: q*("~ mu[1]~")= N("~ .(round(m_iter[i,1],2))~","~
                                                   .(round(sqrt(s2_iter[i,1]),2))~ ")         q*("~ mu[2]~")= N("~ .(round(m_iter[i,2],2))~","~
                                                   .(round(sqrt(s2_iter[i,2]),2))~ ")         q*("~ mu[3]~")= N("~ .(round(m_iter[i,3],2))~","~
                                                   .(round(sqrt(s2_iter[i,3]),2))~ ")"))) +
                  theme(plot.caption = element_text(face = "bold", hjust = 0.5,
                                                    size = rel(1.5)))
            
      })
      
      output$process <- renderPlot({
            combo <- output_combo()
            elbo_iter <- combo$elbo_iter
            psi_iter <- combo$psi_iter
            s2_iter <- combo$s2_iter
            m_iter <- combo$m_iter
            psi <- combo$psi_iter
            i <- combo$i
            
            
            truth <- output_dfsim()
            df_sim <- truth$df_sim
            mu <- truth$mu
            
            # visualize Process
            toplot <- c(1,2,ifelse(i %% 2 == 0, (i/2), (i/2+0.5)), i)
            
            
            # Initial Values
            # generate prediction based on inits
            df_sim_1 <- df_sim %>% 
                  mutate(data = round(data, 5),
                         predicted = as.factor(apply(psi[,,1],1,which.max)))
            # order rows of df_sim by data values
            idx <- order(df_sim_1$data)
            df_sim2_1 <- df_sim_1[idx,]
            binwidth <- diff(range(df_sim_1$data))/35
            
            xlim1 <- c(min(c(m_iter[1,]-1.5*sqrt(s2_iter[1,]),df_sim$data)),
                       max(c(m_iter[1,]+1.5*sqrt(s2_iter[1,]),df_sim$data)))
            plot1 <- df_sim2_1 %>% ggplot(aes(x=data, fill = predicted)) + 
                  geom_dotplot(dotsize = 0.5, binwidth = binwidth,
                               aes(color =  cluster, stroke = 3, alpha = 0.5)) + 
                  theme_bw()+ scale_color_viridis(discrete = TRUE, option = "E") +
                  scale_fill_viridis(discrete = TRUE) + xlim(xlim1) +
                  scale_y_continuous(NULL, breaks = NULL) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[toplot[1],1], sd = sqrt(s2_iter[toplot[1],1])), 
                                colour = 'grey', size = 1, alpha = 0.7, n= 2000) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[toplot[1],2], sd = sqrt(s2_iter[toplot[1],2])), 
                                colour = 'grey', size = 1, alpha = 0.7, n= 2000) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[toplot[1],3], sd = sqrt(s2_iter[toplot[1],3])), 
                                colour = 'grey', size = 1, alpha = 0.7, n= 2000) +
                  geom_vline(xintercept = mu, color = viridis(3, option = "E"), 
                             # linetype = 'dotted',
                             size = 1) +
                  ggtitle("Initial Values") +
                  theme(legend.position = "none") + xlab("") + ylab("") +
                  labs(caption = paste("Assignment Error Rate:", 
                                       100*mean(df_sim_1$cluster !=  as.numeric(df_sim_1$predicted)),
                                       "%")) + 
                  theme(plot.caption = element_text(face = "bold", hjust = 0.5))
            
            
            # first Iter
            # generate prediction based on first iter
            df_sim_2 <- df_sim %>% 
                  mutate(data = round(data, 5),
                         predicted = as.factor(apply(psi[,,2],1,which.max)))
            # order rows of df_sim by data values
            idx <- order(df_sim_2$data)
            df_sim2_2 <- df_sim_2[idx,]
            
            xlim2 <- c(min(c(m_iter[2,]-1.5*sqrt(s2_iter[2,]),df_sim$data)),
                       max(c(m_iter[2,]+1.5*sqrt(s2_iter[2,]),df_sim$data)))
            plot2 <- df_sim2_2 %>% ggplot(aes(x=data, fill = predicted)) + 
                  geom_dotplot(dotsize = 0.5, binwidth = binwidth,
                               aes(color = cluster, stroke = 3, alpha = 0.5)) + 
                  theme_bw() + scale_color_viridis(discrete = TRUE, option = "E") +
                  scale_fill_viridis(discrete = TRUE) +
                  scale_y_continuous(NULL, breaks = NULL) + xlim(xlim2) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[toplot[2],1], sd = sqrt(s2_iter[toplot[2],1])), 
                                colour = 'grey', size = 1, alpha = 0.7, n= 2000) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[toplot[2],2], sd = sqrt(s2_iter[toplot[2],2])), 
                                colour = 'grey', size = 1, alpha = 0.7, n= 2000) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[toplot[2],3], sd = sqrt(s2_iter[toplot[2],3])), 
                                colour = 'grey', size = 1, alpha = 0.7, n= 2000) +
                  geom_vline(xintercept = mu, color = viridis(3, option = "E"), 
                             # linetype = 'dotted',
                             size = 1) +
                  ggtitle("First Iteration") +
                  theme(legend.position = "none") + xlab("") + ylab("") +
                  labs(caption = paste("Assignment Error Rate:", 
                                       100*mean(df_sim_2$cluster !=  as.numeric(df_sim_2$predicted)),
                                       "%")) + 
                  theme(plot.caption = element_text(face = "bold", hjust = 0.5))
            
            # halfway Iter
            # generate prediction based on halfway iter
            df_sim_3 <- df_sim %>% 
                  mutate(data = round(data, 5),
                         predicted = as.factor(apply(psi[,,toplot[3]],1,which.max)))
            # order rows of df_sim by data values
            idx <- order(df_sim_3$data)
            df_sim2_3 <- df_sim_3[idx,]
            
            xlim3 <- c(min(c(m_iter[toplot[3],]-1.5*sqrt(s2_iter[toplot[3],]),df_sim$data)),
                       max(c(m_iter[toplot[3],]+1.5*sqrt(s2_iter[toplot[3],]),df_sim$data)))
            plot3 <- df_sim2_3 %>% ggplot(aes(x=data, fill = predicted)) + 
                  geom_dotplot(dotsize = 0.5, binwidth = binwidth,
                               aes(color = cluster, stroke = 3, alpha = 0.5)) + 
                  theme_bw() + scale_color_viridis(discrete = TRUE, option = "E") +
                  scale_fill_viridis(discrete = TRUE) +
                  scale_y_continuous(NULL, breaks = NULL) + xlim(xlim3) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[toplot[3],1], 
                                                         sd = sqrt(s2_iter[toplot[3],1])), 
                                colour = 'grey', size = 1, alpha = 0.7, n= 2000) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[toplot[3],2], 
                                                         sd = sqrt(s2_iter[toplot[3],2])), 
                                colour = 'grey', size = 1, alpha = 0.7, n= 2000) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[toplot[3],3], 
                                                         sd = sqrt(s2_iter[toplot[3],3])), 
                                colour = 'grey', size = 1, alpha = 0.7, n= 2000) +
                  geom_vline(xintercept = mu, color = viridis(3, option = "E"), 
                             # linetype = 'dotted',
                             size = 1) +
                  ggtitle(paste("After ",toplot[3]," Iterations", sep = ""))+
                  theme(legend.position = "none") + xlab("") + ylab("") +
                  labs(caption = paste("Assignment Error Rate:", 
                                       100*mean(df_sim_3$cluster !=  as.numeric(df_sim_3$predicted)),
                                       "%")) + 
                  theme(plot.caption = element_text(face = "bold", hjust = 0.5))
            
            # Last Iter
            # generate prediction based on last iter
            df_sim_4 <- df_sim %>% 
                  mutate(data = round(data, 5),
                         predicted = as.factor(apply(psi[,,i],1,which.max)))
            # order rows of df_sim by data values
            idx <- order(df_sim_4$data)
            df_sim2_4 <- df_sim_4[idx,]
            
            plot4 <- df_sim2_4 %>% ggplot(aes(x=data, fill = predicted)) + 
                  geom_dotplot(dotsize = 0.5, binwidth = binwidth,
                               aes(color = cluster, stroke = 3, alpha = 0.5)) + 
                  theme_bw() + scale_color_viridis(discrete = TRUE, option = "E")+
                  scale_fill_viridis(discrete = TRUE) +
                  scale_y_continuous(NULL, breaks = NULL) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[toplot[4],1], sd = sqrt(s2_iter[toplot[4],1])), 
                                colour = 'grey', size = 1, alpha = 0.7, n= 2000) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[toplot[4],2], sd = sqrt(s2_iter[toplot[4],2])), 
                                colour = 'grey', size = 1, alpha = 0.7, n= 2000) +
                  stat_function(fun = dnorm, args = list(mean =m_iter[toplot[4],3], sd = sqrt(s2_iter[toplot[4],3])), 
                                colour = 'grey', size = 1, alpha = 0.7, n= 2000) +
                  geom_vline(xintercept = mu, color = viridis(3, option = "E"), 
                             # linetype = 'dotted',
                             size = 1) +
                  ggtitle(paste("After ", toplot[4]," Iterations", sep = ""))+
                  theme(legend.position = "none") + xlab("") + ylab("") +
                  labs(caption = paste("Assignment Error Rate:", 
                                       100*mean(df_sim_4$cluster !=  as.numeric(df_sim_4$predicted)),
                                       "%")) + 
                  theme(plot.caption = element_text(face = "bold", hjust = 0.5))
            
            
            grid.arrange(plot1, plot2, plot3, plot4, 
                         nrow = 2, ncol = 2)
            
      })
      output$elbo <- renderPlot({
            combo <- output_combo()
            elbo_iter <- combo$elbo_iter
            i <- combo$i
            elbo_iter_plot <- as.data.frame(elbo_iter[1:i]) %>% mutate(iter = 1:i)
            elbo_iter_plot %>% ggplot(aes(x = iter, y = elbo_iter[1:i])) + 
                  geom_point() + geom_smooth(se = FALSE) + theme_bw() + xlab("Iteration")  +
                  scale_x_continuous(breaks = seq(1:i)) +
                  ylab("ELBO")
      })
      
      output$table <- DT::renderDataTable({
            truth <- output_dfsim()
            df_sim <- truth$df_sim
            combo <- output_combo()
            psi <- combo$psi_iter
            i <- combo$i
            df_sim <- df_sim %>%
                  mutate(data = round(data, 3),
                         cluster = as.numeric(cluster),
                         predicted = apply(psi[,,i],1,which.max),
                         match = ifelse(predicted == cluster, "yes", "no"),
                         psi1 = round(psi[,1,i],2),
                         psi2 = round(psi[,2,i],2),
                         psi3 = round(psi[,3,i],2))
            colnames(df_sim) <- c("Simulated Values",
                                  "True Assignment",
                                  "Predicted Assignment",
                                  "Assignment Correctly Predicted?",
                                  "Predicted Assignment Probability for Group1",
                                  "Predicted Assignment Probability for Group2",
                                  "Predicted Assignment Probability for Group3")
            df_sim
      })
      
}

# Run the application 
shinyApp(ui = ui, server = server)


