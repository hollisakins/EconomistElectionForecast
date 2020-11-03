
#' Description
#' This file runs a live election-night forecast based on The Economist's pre-election forecasting model
#' available at projects.economist.com/us-2020-forecast/president.
#' It is resampling model based on https://pkremp.github.io/update_prob.html.
#' This script does not input any real election results! You will have to enter your picks/constraints manually 
#' (scroll to the bottom of the script).
#' 
#' Modified by Hollis Akins (@hollisakins) to run in a shiny app. 
#' Not the fastest or most efficient code (and not commented that well) but it works just fine
#' 

library(shiny)
library(tidyverse)
library(mvtnorm)
library(politicaldata)
library(ggplot2)
library(usmap)
library(ggpubr)

# functions first ---------------------------------------------------------
logit <- function(x) log(x/(1-x))

inv_logit <- function(x) 1/(1 + exp(-x))


draw_samples <- function(biden_states = NULL, trump_states = NULL, states = NULL, 
                         upper_biden = NULL, lower_biden = NULL, print_acceptance = FALSE, target_nsim = 1000){
    sim <- matrix(NA, nr = 1, nc = length(mu))
    n <- 0
    while(nrow(sim) < target_nsim){
        # randomly sample from the posterior distribution and reject when constraints are not met
        n <- n + 1
        proposals <- inv_logit(rmvnorm(1e5, mu, Sigma, method = "svd")) # "DC" is pretty much uncorrelated
        colnames(proposals) <- names(mu)
        if (!is.null(biden_states)) proposals[which(proposals[,biden_states] < .5)] <- NA
        if (!is.null(  trump_states)) proposals[which(proposals[,  trump_states] > .5)] <- NA
        if (!is.null(        states)){
            for (s in states){
                proposals[which(proposals[, s] > upper_biden[s] | 
                                    proposals[, s] < lower_biden[s])] <- NA
            }
        }
        reject <- apply(proposals, 1, function(x) any(is.na(x)))
        sim <- rbind(sim, proposals[!reject,])
        if (nrow(sim) < target_nsim & nrow(sim)/(nrow(proposals)*n) < 1-99.99/100){
            stop(paste("rmvnorm() is working hard... but more than 99.99% of the samples are rejected; you should relax some contraints.", sep = ""))
        }
    }
    return(list("matrix" = sim[-1,], "acceptance_rate" = nrow(sim)/(nrow(proposals)*n)))
}


update_prob <- function(biden_states = NULL, trump_states = NULL, biden_scores_list = NULL, target_nsim = 1000, show_all_states = TRUE, main_plots = TRUE, side_plots = FALSE){
    
    states <- names(biden_scores_list)
    lower_biden <- sapply(biden_scores_list, function(x) x[1]/100)
    upper_biden <- sapply(biden_scores_list, function(x) x[2]/100)
    sim <- draw_samples(biden_states = biden_states, trump_states = trump_states, states = states, 
                        upper_biden = upper_biden, lower_biden = lower_biden, 
                        target_nsim = target_nsim)
    ev_dist <- (sim[["matrix"]] > .5) %*% ev
    
    state_win <- colMeans(sim[["matrix"]] > .5)
    p <- mean(ev_dist >= 270)
    sd <- sqrt(p*(1-p)/length(ev_dist))
    
    ev_df <- as.data.frame(ev)
    ev_df$state <- row.names(ev_df)
    
    state_win_df <- as.data.frame(state_win)
    state_win_df$state <- row.names(state_win_df)
    
    df <- dplyr::full_join(ev_df,state_win_df,by="state")
    df$dummy <- "dummy"
    df$dummy <- as.factor(df$dummy)
    
    df$color <- as.factor(ifelse(df$state %in% biden_states | df$state %in% trump_states, 
                       yes = "white", no = "black"))
    
    df <- df[order(df$state_win),]
    df$state <- ordered(df$state)
    df$color <- ordered(df$color)
    
    # if (show_all_states){
    #     cat("Pr(biden wins) by state, in %:\n")
    #     print(t(round(100*state_win,1)))
    #     cat("--------\n")
    # }
    
    #cat(paste("Pr(biden wins the electoral college) = ", round(100*p,1), "%\n[nsim = ", length(ev_dist), "; se = ", round(sd*100,1), "%]", sep = ""))
    #if (show_all_states) cat("\n--------\n")
    
    if (length(biden_states)==0 & length(trump_states)==0) {
        tit <- "Economist forecast"
    }
    if (!(length(biden_states)==0) & length(trump_states)==0) {
        tit <- paste("Biden wins",biden_states)
    }
    if (!(length(trump_states)==0) & is.null(biden_states)) {
        tit <- paste("Trump wins",trump_states)
    }
    if (!(length(trump_states)==0) & !(length(biden_states)==0)) {
        tit <- paste("Biden wins",biden_states,", Trump wins",trump_states)
    }
    
    
    
    
    
    if (main_plots){
        
       g1 <- plot_usmap(data=state_win_df, values="state_win", labels=T) + 
          scale_fill_gradient2(midpoint=0.5,name ="Biden Win Prob") +
          theme(legend.position = "right", 
                plot.margin=grid::unit(c(0,0,0,0), "mm"),
                text = element_text(size=20)) + 
          labs(title=tit) 
    
        g2 <- ggplot(data = df, aes(x=ev, y=dummy, fill=state_win))+#, color=color)) + 
            geom_bar(stat="identity", color="black") + 
            geom_vline(xintercept=270, linetype="dotted", color="black", size=1) + 
            xlab("270 to win") +
            scale_fill_gradient2(midpoint=0.5, aesthetics="fill") +
            #scale_color_manual(values=c("black","white"), aesthetics = "color") +
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.y=element_blank(),legend.position="none",
                      panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(),
                  text = element_text(size=15))
            
        gfinal <- ggarrange(g1,g2,ncol = 1, nrow = 2, heights=c(6,1))
        return(gfinal)
    }
    
    if (side_plots){
        text <- paste0("Biden Win Probability = ", round(100*p,1), "%")
        
        gfinal <- ggplot(data=data.frame(ev=ev_dist)) +
            geom_histogram(aes(x=ev), bins=538) + 
            geom_density(aes(x=ev, y = ..count..), color="red", size=1.5) + 
            theme(plot.background = element_rect(fill = "transparent",colour = NA),
                  text = element_text(size=15)) + 
            labs(title=text) + 
            xlab("Electoral College Total") + 
            ylab("Count") + xlim(0, 538) + 
            geom_vline(xintercept=270, linetype="dotted", color="black", size=1)
        
        return(gfinal)
    }
   
}


# read data ---------------------------------------------------------------
# get simulations
sim_forecast <- read_csv('https://cdn.economistdatateam.com/us-2020-forecast/data/president/electoral_college_simulations.csv') 

# check initial parameters
nrow(sim_forecast)
mean(sim_forecast$dem_ev > 269)

# select relevant columns and make the data frae into a matrix
sim_forecast <- sim_forecast %>% select(4:ncol(.))
sim_forecast <- list(sim_forecast,sim_forecast,sim_forecast,sim_forecast,sim_forecast) %>% bind_rows  %>% as.matrix


# this bit is really hacky
# it make the simulations a little less correlated by add a bunch of random noise
# this helps our live model react to really implausible events that could happen on election night
# but will have the result of pushing the initial pre-results forecasts back toward 50-50
sim_forecast <- sim_forecast +
    rnorm(nrow(sim_forecast), 0, 0.01) +  # national error component
    replicate(ncol(sim_forecast),rnorm(nrow(sim_forecast),0,0.02)) # state

sim_forecast <- ifelse(sim_forecast <= 0, 0.0001, sim_forecast)
sim_forecast <- ifelse(sim_forecast >= 1, 0.99999, sim_forecast)

sim_forecast <- as_tibble(sim_forecast)

# now, get electoral votes in each state 
# and make sure they're in the right order
#ev_state <- read_csv('data/state_evs.csv')$ev
#names(ev_state) <- read_csv('data/state_evs.csv')$state
ev_state <- read_csv('https://raw.githubusercontent.com/TheEconomist/us-potus-model/master/data/2012.csv')$ev
names(ev_state) <- read_csv('https://raw.githubusercontent.com/TheEconomist/us-potus-model/master/data/2012.csv')$state
ev <- ev_state[colnames(sim_forecast)] 

# check that the EVs and state forecasts add up to the right amounts
sim_evs <- apply(sim_forecast,
                 1,
                 function(x){
                     sum((x > 0.5 )* ev)
                 })

#hist(sim_evs,breaks=538) 
#median(sim_evs)
#mean(sim_evs)
#mean(sim_evs > 269)
#enframe(prop.table(table(sim_evs)),'dem_ev','pct') %>% arrange(desc(pct))


# adding ME1 and ME2, NE1 NE2 to sim_forecast matrix and ev vector
# we do this by adding the average district-level two-party dem presidential vote, relative
# to the state-level dem two-party vote, back to to our state-level forecast

# first, split up EVs
ev["ME"] <- 2
ev["NE"] <- 2
ev <- c(ev, "ME1" = 1, "ME2" = 1, "NE1" = 1, "NE2" = 1, "NE3" = 1)
sum(ev)

# create simulations for ME and NE districts
me_ne_leans <- politicaldata::pres_results_by_cd %>% filter(year >= 2012, state_abb %in% c("ME","NE")) %>%
    select(-other) %>% 
    rename(state = state_abb) %>%
    group_by(year,state) %>%
    mutate(sum_pct = dem + rep,
           total_votes = total_votes * sum_pct,
           dem = dem /sum_pct,
           rep = rep/sum_pct) %>%
    mutate(dem_vote_state = sum(dem * total_votes) / sum(total_votes)) %>%
    mutate(dem_cd_lean = dem - dem_vote_state) %>%
    group_by(state,district) %>%
    summarise(dem_cd_lean = weighted.mean(dem_cd_lean, c(0.3,0.7)))

# bind new simulation columns for the congressional districts, based on the above
sim_forecast <- bind_cols(
    sim_forecast, 
    tibble(
        ME1 = sim_forecast %>% pull(ME) +
            rnorm(nrow(sim_forecast), 
                  me_ne_leans[me_ne_leans$state == "ME" & me_ne_leans$district == 1,]$dem_cd_lean, 
                  .0075),
        ME2 = sim_forecast %>% pull(ME) +
            rnorm(nrow(sim_forecast), 
                  me_ne_leans[me_ne_leans$state == "ME" & me_ne_leans$district == 2,]$dem_cd_lean, 
                  .0075),
        
        NE1 = sim_forecast %>% pull(NE) +
            rnorm(nrow(sim_forecast), 
                  me_ne_leans[me_ne_leans$state == "NE" & me_ne_leans$district == 1,]$dem_cd_lean, 
                  .0075),
        NE2 = sim_forecast %>% pull(NE) +
            rnorm(nrow(sim_forecast), 
                  me_ne_leans[me_ne_leans$state == "NE" & me_ne_leans$district == 2,]$dem_cd_lean, 
                  .0075),
        NE3 = sim_forecast %>% pull(NE) +
            rnorm(nrow(sim_forecast), 
                  me_ne_leans[me_ne_leans$state == "NE" & me_ne_leans$district == 3,]$dem_cd_lean, 
                  .0075) 
    )
)

# sim_forecast


sim_evs <- apply(sim_forecast,
                 1,
                 function(x){
                     sum((x > 0.5 )* ev)
                 })

# colMeans(sim_forecast > 0.5)
# 
# hist(sim_evs,breaks=538)
# median(sim_evs)
# mean(sim_evs)
# mean(sim_evs > 269)
# enframe(prop.table(table(sim_evs)),'dem_ev','pct') %>% arrange(desc(pct))


# final data wrangling -- this stuff getes passed into the update_prob function
# mainly we just want to make sure everything is in the right order with the right names
ev <- ev[colnames(sim_forecast)]

Sigma <- cov(logit(sim_forecast)) 
Sigma

mu <- colMeans(logit(sim_forecast))
names(mu) <- colnames(sim_forecast)



# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("2020 Election Forecasting, via the Economist"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("bs", strong("Biden wins:"), choices = state.abb, selected = NULL, multiple = T),
            selectInput("ts", strong("Trump wins:"), choices = state.abb, selected = NULL, multiple = T),
            plotOutput("sideplots", width="100%", height="400px")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("map", width="100%", height="800px")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    
    
    output$map <- renderPlot({
        
        g <- update_prob(biden_states = input$bs, 
                         trump_states = input$ts,
                         biden_scores_list = NULL)
        print(g)})
    
    output$sideplots <- renderPlot({
        g <- update_prob(biden_states = input$bs,
                         trump_states = input$ts,
                         biden_scores_list = NULL, 
                         side_plots = T, main_plots = F)
        print(g)}, bg="transparent")
}

# Run the application 
shinyApp(ui = ui, server = server)
