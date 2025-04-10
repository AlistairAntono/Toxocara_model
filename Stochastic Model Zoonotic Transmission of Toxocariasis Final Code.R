library(dplyr)
library(tidyr)
library(optimx)

data.frame(Input_parameters)

#########################################
### MAIN MODEL FUNCTION #################
#########################################

run_country_model <- function(country_name, Input_parameters, Input_demographics, transmission_params) {
  # Progress message
  print(paste("Processing", country_name))
  
  # Extract parameters for country
  country_params <- Input_parameters[Input_parameters$Country == country_name, ]
  country_demos <- data.frame(
    Age = Input_demographics$Age_group,
    Population = as.numeric(unlist(Input_demographics[, country_name]))
  )
  
  # Get region for this country and corresponding transmission parameters
  country_region <- country_params$Region
  
  # Create transmission parameters list using regional values
  transmission_params <- list(
    beta1 = transmission_params_long$beta1[transmission_params_long$Parameter_Region == country_region],
    beta2 = transmission_params_long$beta2[transmission_params_long$Parameter_Region == country_region],
    beta1_ci_lower = transmission_params_long$beta1_ci_lower[transmission_params_long$Parameter_Region == country_region],
    beta1_ci_upper = transmission_params_long$beta1_ci_upper[transmission_params_long$Parameter_Region == country_region],
    beta2_ci_lower = transmission_params_long$beta2_ci_lower[transmission_params_long$Parameter_Region == country_region],
    beta2_ci_upper = transmission_params_long$beta2_ci_upper[transmission_params_long$Parameter_Region == country_region]
  )
  
  # Extract all other parameters
  animal1_pop_2021 <- country_params$animal1_pop_2021
  animal1_pop_2022 <- country_params$animal1_pop_2022
  animal1_pop_2023 <- country_params$animal1_pop_2023
  animal2_pop_2021 <- country_params$animal2_pop_2021
  animal2_pop_2022 <- country_params$animal2_pop_2022
  animal2_pop_2023 <- country_params$animal2_pop_2023
  initial_population_human <- country_params$initial_population_human
  human_growth_rate_yearly <- country_params$human_growth_rate_yearly
  growth_rate_human <- human_growth_rate_yearly / 12
  animal1_prevalence <- country_params$animal1_prevalence
  animal2_prevalence <- country_params$animal2_prevalence
  initial_recovered_proportion <- country_params$initial_recovered_proportion
  estimated_infections <- country_params$estimated_infections
  covert_cost <- country_params$covert_cost
  ocular_cost <- country_params$ocular_cost
  visceral_cost <- country_params$visceral_cost
  
  # Calculate growth rates
  calculate_growth_rate <- function(pop_2021, pop_2022, pop_2023, max_reasonable_rate = 0.05) {
    growth_rate_2021_2022 <- (pop_2022 - pop_2021) / pop_2021
    growth_rate_2022_2023 <- (pop_2023 - pop_2022) / pop_2022
    selected_rate <- min(abs(growth_rate_2021_2022), abs(growth_rate_2022_2023))
    final_rate <- min(selected_rate, max_reasonable_rate)
    return(final_rate)
  }
  
  growth_rate_animal1 <- calculate_growth_rate(animal1_pop_2021, animal1_pop_2022, animal1_pop_2023)
  growth_rate_animal2 <- calculate_growth_rate(animal2_pop_2021, animal2_pop_2022, animal2_pop_2023)
  growth_rate_animal1_monthly <- growth_rate_animal1 / 12
  growth_rate_animal2_monthly <- growth_rate_animal2 / 12
  
  # Calculate transmission parameters and CIs
  percentage_cases_age_group <- c(14, 10, 16, 12, 16, 12, 9, 8, 3)
  populations <- country_demos$Population[1:9]
  total_population <- country_demos$Population[10]
  percentage_population_age_group <- (populations / total_population) * 100
  scaling_factors <- percentage_cases_age_group / percentage_population_age_group
  
  
  # Save initial parameters including CIs
  param_df <- data.frame(
    Parameter = c("beta1", "beta2", 
                  "beta1_ci_lower", "beta1_ci_upper",
                  "beta2_ci_lower", "beta2_ci_upper",
                  "growth_rate_animal1_yearly", "growth_rate_animal2_yearly",
                  "growth_rate_animal1_monthly", "growth_rate_animal2_monthly",
                  paste0("scaling_factor_age_group_", 1:9)),
    Value = c(transmission_params$beta1, transmission_params$beta2,
              transmission_params$beta1_ci_lower, transmission_params$beta1_ci_upper,
              transmission_params$beta2_ci_lower, transmission_params$beta2_ci_upper,
              growth_rate_animal1, growth_rate_animal2,
              growth_rate_animal1_monthly, growth_rate_animal2_monthly,
              scaling_factors)
  )
  write.csv(param_df, paste0(country_name, "_parameters.csv"), row.names = FALSE)
  
  rpert <- function(n, min, max, mode, lambda = 4) {
    # alpha1 and alpha2 parameters for the beta distribution
    alpha1 <- 1 + lambda * (mode - min)/(max - min)
    alpha2 <- 1 + lambda * (max - mode)/(max - min)
    
    # Generate samples from beta distribution
    x <- rbeta(n, alpha1, alpha2)
    
    # Transform to desired range
    min + (max - min) * x
  }
  
  # Initialize storage for all iterations
  all_iterations_results <- list()
  
  #########################################
  ### ITERATION LOOP #####################
  #########################################
  
  for(iteration in 1:1000) {
    # Set seed for this iteration
    set.seed(iteration * 1000)
    
    # Fixed parameters
    gamma <- 0.278 / 12
    reversion_rate <- 1 / (12 * 1)
    treatment_effectiveness1 <- 0.0 / 12
    treatment_effectiveness2 <- 0.0 / 12
    
    # Initialize populations dataframe with beta columns
    populations <- data.frame(
      Month = 1:120,
      Susceptible_Humans = numeric(120),
      Infected_Humans = numeric(120),
      Recovered_Humans = numeric(120),
      Animals1 = numeric(120),
      Animals2 = numeric(120),
      Total_New_Exposures = numeric(120),  
      Symptomatic_Cases = numeric(120),   
      Severe_Cases = numeric(120),
      Covert_Cases = numeric(120),
      Ocular_Cases = numeric(120),
      Visceral_Cases = numeric(120),
      Monthly_Cost_Covert = numeric(120),
      Monthly_Cost_Ocular = numeric(120),
      Monthly_Cost_Visceral = numeric(120),
      Cumulative_Total_Cost = numeric(120),
      beta1_monthly = numeric(120),
      beta2_monthly = numeric(120),
      Cumulative_Exposures = numeric(120) 
    )
    
    # Set initial populations
    initial_recovered <- initial_population_human * initial_recovered_proportion
    populations[1, c("Susceptible_Humans", "Recovered_Humans")] <- 
      c(initial_population_human - initial_recovered, initial_recovered)
    populations[1, c("Animals1", "Animals2")] <- 
      c(animal1_pop_2023, animal2_pop_2023)
    
    # Initialize storage for cumulative counts
    cumulative_exposures <- 0
    populations$Cumulative_Exposures <- numeric(120)
    
    #########################################
    ### MONTHLY SIMULATION LOOP #############
    #########################################
    
    for (i in 2:120) {
      # Sample this month's beta values from within CIs with validation
      beta1_ci_lower <- ifelse(is.infinite(transmission_params$beta1_ci_lower), 
                               transmission_params$beta1 * 0.8, 
                               transmission_params$beta1_ci_lower)
      beta1_ci_upper <- ifelse(is.infinite(transmission_params$beta1_ci_upper), 
                               transmission_params$beta1 * 1.2, 
                               transmission_params$beta1_ci_upper)
      beta2_ci_lower <- ifelse(is.infinite(transmission_params$beta2_ci_lower), 
                               transmission_params$beta2 * 0.8, 
                               transmission_params$beta2_ci_lower)
      beta2_ci_upper <- ifelse(is.infinite(transmission_params$beta2_ci_upper), 
                               transmission_params$beta2 * 1.2, 
                               transmission_params$beta2_ci_upper)
      
      beta1_month <- rpert(1, beta1_ci_lower, beta1_ci_upper, transmission_params$beta1)
      beta2_month <- rpert(1, beta2_ci_lower, beta2_ci_upper, transmission_params$beta2)
      
      # Store sampled beta values
      populations$beta1_monthly[i] <- beta1_month
      populations$beta2_monthly[i] <- beta2_month
      
      # Calculate exposures using this month's betas
      total_new_exposures_age_group <- sapply(1:length(scaling_factors), function(j) {
        beta_age_group1 <- beta1_month * scaling_factors[j]
        beta_age_group2 <- beta2_month * scaling_factors[j]
        new_exposures1 <- beta_age_group1 * animal1_prevalence * 
          populations[i - 1, "Susceptible_Humans"] * 
          populations[i - 1, "Animals1"] / initial_population_human
        new_exposures2 <- beta_age_group2 * animal2_prevalence * 
          populations[i - 1, "Susceptible_Humans"] * 
          populations[i - 1, "Animals2"] / initial_population_human
        return(new_exposures1 + new_exposures2)
      })
      
      total_new_exposures <- sum(total_new_exposures_age_group)
      
      # Calculate symptomatic cases with stochastic proportion using PERT
      symptomatic_proportion <- rpert(1, min = 0.20, max = 0.40, mode = 0.30)
      populations$Symptomatic_Cases[i] <- total_new_exposures * symptomatic_proportion
      
      # Update populations
      populations[i, "Total_New_Exposures"] <- total_new_exposures
      populations$Cumulative_Exposures[i] <- cumulative_exposures + total_new_exposures
      cumulative_exposures <- populations$Cumulative_Exposures[i]
      
      # Update human compartments
      new_exposed <- total_new_exposures
      # Update susceptible (reduced by exposures)
      populations[i, "Susceptible_Humans"] <- populations[i-1, "Susceptible_Humans"] - new_exposed
      
      # Update infected (increased by exposures, reduced by recovery)
      populations[i, "Infected_Humans"] <- populations[i-1, "Infected_Humans"] + new_exposed - 
        (populations[i-1, "Infected_Humans"] * gamma)
      
      # Update recovered (increased by recovery from infected)
      populations[i, "Recovered_Humans"] <- populations[i-1, "Recovered_Humans"] + 
        (populations[i-1, "Infected_Humans"] * gamma)
      
      # Update animal populations with growth
      populations[i, "Animals1"] <- populations[i-1, "Animals1"] * (1 + growth_rate_animal1_monthly)
      populations[i, "Animals2"] <- populations[i-1, "Animals2"] * (1 + growth_rate_animal2_monthly)
      
      # Human Population Growth (after all other population updates)
      populations[i, c("Susceptible_Humans", "Infected_Humans", "Recovered_Humans")] <- 
        populations[i, c("Susceptible_Humans", "Infected_Humans", "Recovered_Humans")] * 
        (1 + growth_rate_human)
      
      # Calculate case numbers based on symptomatic cases
      severe_cases_percentage <- (20 / 12) / populations[2, "Symptomatic_Cases"]
      populations$Severe_Cases[i] <- populations$Symptomatic_Cases[i] * severe_cases_percentage
      populations$Covert_Cases[i] <- populations$Symptomatic_Cases[i] - populations$Severe_Cases[i]
      populations$Ocular_Cases[i] <- populations$Severe_Cases[i] * 0.21
      populations$Visceral_Cases[i] <- populations$Severe_Cases[i] * 0.79
      
      # Calculate costs
      populations$Monthly_Cost_Covert[i] <- populations$Covert_Cases[i] * 0.3 * covert_cost 
      populations$Monthly_Cost_Ocular[i] <- populations$Ocular_Cases[i] * ocular_cost
      populations$Monthly_Cost_Visceral[i] <- populations$Visceral_Cases[i] * visceral_cost
      
      # Calculate cumulative cost
      populations$Cumulative_Total_Cost[i] <- sum(populations$Monthly_Cost_Covert[2:i]) +
        sum(populations$Monthly_Cost_Ocular[2:i]) +
        sum(populations$Monthly_Cost_Visceral[2:i])
    }
    
    # Add iteration number to populations dataframe
    populations$Iteration <- iteration
    
    # Save full results for this iteration
    if(iteration == 1) {
      write.csv(populations, paste0(country_name, "_populations.csv"), row.names = FALSE)
    } else {
      write.table(populations, paste0(country_name, "_populations.csv"), 
                  append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
    }
    
    # Store iteration results
    iteration_result <- data.frame(
      Country = country_name,
      Iteration = iteration,
      Cumulative_Exposures = tail(populations$Cumulative_Exposures, 1),
      Total_Symptomatic = sum(populations$Symptomatic_Cases),
      Severe_Cases = sum(populations$Severe_Cases),
      Covert_Cases = sum(populations$Covert_Cases),
      Ocular_Cases = sum(populations$Ocular_Cases),
      Visceral_Cases = sum(populations$Visceral_Cases),
      Cumulative_Total_Cost = tail(populations$Cumulative_Total_Cost, 1),
      Mean_Beta1 = mean(populations$beta1_monthly, na.rm = TRUE),
      Mean_Beta2 = mean(populations$beta2_monthly, na.rm = TRUE)
    )
    
    all_iterations_results[[iteration]] <- iteration_result
  }  # End of iteration loop
  
  # Combine all iterations
  country_results <- do.call(rbind, all_iterations_results)
  
  # Calculate summary statistics across iterations
  summary_stats <- data.frame(
    Country = country_name,
    Iteration = "Summary",
    Cumulative_Exposures_Mean = mean(country_results$Cumulative_Exposures),
    Cumulative_Exposures_SD = sd(country_results$Cumulative_Exposures),
    Total_Symptomatic_Mean = mean(country_results$Total_Symptomatic),
    Total_Symptomatic_SD = sd(country_results$Total_Symptomatic),
    Severe_Cases_Mean = mean(country_results$Severe_Cases),
    Severe_Cases_SD = sd(country_results$Severe_Cases),
    Covert_Cases_Mean = mean(country_results$Covert_Cases),
    Covert_Cases_SD = sd(country_results$Covert_Cases),
    Ocular_Cases_Mean = mean(country_results$Ocular_Cases),
    Ocular_Cases_SD = sd(country_results$Ocular_Cases),
    Visceral_Cases_Mean = mean(country_results$Visceral_Cases),
    Visceral_Cases_SD = sd(country_results$Visceral_Cases),
    Cumulative_Total_Cost_Mean = mean(country_results$Cumulative_Total_Cost),
    Cumulative_Total_Cost_SD = sd(country_results$Cumulative_Total_Cost),
    Beta1_Mean = mean(country_results$Mean_Beta1),
    Beta1_SD = sd(country_results$Mean_Beta1),
    Beta2_Mean = mean(country_results$Mean_Beta2),
    Beta2_SD = sd(country_results$Mean_Beta2)
  )
  
  print(paste("Complete", country_name))
  
  # Return both iteration results and summary statistics
  return(list(iterations = country_results, summary = summary_stats))
}
#########################################
### MAIN PROCESSING LOOP ################
#########################################


# Clear and improved main processing section
print("Starting main processing loop")
print("Countries to process:")
print(Input_parameters$Country)

# Reshape transmission parameters before processing
transmission_params_long <- transmission_params %>%
  pivot_longer(
    cols = c(Europe, `South America`, Africa, `North America`, Asia),
    names_to = "Parameter_Region",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Parameter,
    values_from = Value
  )

# Process all countries with error handling
all_results <- lapply(Input_parameters$Country, function(country) {
  print(paste("Attempting to process:", country))
  tryCatch({
    result <- run_country_model(country, Input_parameters, Input_demographics, transmission_params_long)
    print(paste("Successfully processed:", country))
    return(result)
  }, error = function(e) {
    print(paste("Error processing", country, ":", e$message))
    return(NULL)
  })
})
# Remove any NULL results from failed processing
all_results <- all_results[!sapply(all_results, is.null)]

# Before combining results, check if we have any valid results
if (length(all_results) > 0) {
  print("Number of successful results:")
  print(length(all_results))
  print("First result structure:")
  print(str(all_results[[1]]))
  
  # Combine all results
  final_results_iterations <- do.call(rbind, lapply(all_results, function(x) x$iterations))
  final_results_summary <- do.call(rbind, lapply(all_results, function(x) x$summary))
  
  # Save both iterations and summary results
  write.csv(final_results_iterations, "Final_results_iterations.csv", row.names = FALSE)
  write.csv(final_results_summary, "Final_results_summary.csv", row.names = FALSE)
  
  print("Results have been saved to CSV files")
} else {
  print("No successful results to process")
}


#########################################
### MULTI PANEL PLOTS ###################
#########################################

create_country_summary_plot <- function(country_name) {
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  
  # Read the populations data
  populations <- read.csv(paste0(country_name, "_populations.csv"))
  
  # Calculate summary statistics by month
  summary_by_month <- populations %>%
    filter(Month >= 2) %>%
    group_by(Month) %>%
    summarise(
      # Total Exposures
      mean_exposures = mean(Total_New_Exposures),
      lower_ci_exposures = quantile(Total_New_Exposures, 0.025),
      upper_ci_exposures = quantile(Total_New_Exposures, 0.975),
      # Symptomatic Cases
      mean_symptomatic = mean(Symptomatic_Cases),
      lower_ci_symptomatic = quantile(Symptomatic_Cases, 0.025),
      upper_ci_symptomatic = quantile(Symptomatic_Cases, 0.975),
      # Cumulative exposures
      mean_cumulative = mean(Cumulative_Exposures),
      lower_ci_cumulative = quantile(Cumulative_Exposures, 0.025),
      upper_ci_cumulative = quantile(Cumulative_Exposures, 0.975),
      # Costs
      mean_cost = mean(Cumulative_Total_Cost),
      lower_ci_cost = quantile(Cumulative_Total_Cost, 0.025),
      upper_ci_cost = quantile(Cumulative_Total_Cost, 0.975),
      # Animal populations and prevalence
      mean_animals1 = mean(Animals1),
      mean_animals2 = mean(Animals2)
    )
  
  # Plot 1: Cumulative Exposures
  p1 <- ggplot(summary_by_month, aes(x = Month)) +
    geom_ribbon(aes(ymin = lower_ci_cumulative, ymax = upper_ci_cumulative),
                fill = "lightblue", alpha = 0.3) +
    geom_line(aes(y = mean_cumulative), color = "darkblue") +
    labs(title = "Cumulative Exposures Over Time",
         y = "Cumulative Exposures") +
    theme_minimal()
  
  # Plot 2: Monthly Exposures and Cases
  p2 <- ggplot(summary_by_month) +
    geom_ribbon(aes(x = Month, ymin = lower_ci_exposures, ymax = upper_ci_exposures),
                fill = "lightblue", alpha = 0.3) +
    geom_line(aes(x = Month, y = mean_exposures), color = "darkblue") +
    geom_ribbon(aes(x = Month, ymin = lower_ci_symptomatic, ymax = upper_ci_symptomatic),
                fill = "pink", alpha = 0.3) +
    geom_line(aes(x = Month, y = mean_symptomatic), color = "darkred") +
    labs(title = "Monthly Exposures and Cases",
         y = "Cases per Month") +
    theme_minimal()
  
  # Plot 3: Animal Populations
  p3 <- ggplot(summary_by_month, aes(x = Month)) +
    geom_line(aes(y = mean_animals1, color = "Animal 1")) +
    geom_line(aes(y = mean_animals2, color = "Animal 2")) +
    scale_color_manual(values = c("Animal 1" = "forestgreen", "Animal 2" = "orange")) +
    labs(title = "Animal Populations Over Time",
         y = "Population Size",
         color = "Species") +
    theme_minimal()
  
  # Plot 4: Cumulative Costs
  p4 <- ggplot(summary_by_month, aes(x = Month)) +
    geom_ribbon(aes(ymin = lower_ci_cost, ymax = upper_ci_cost),
                fill = "lightgreen", alpha = 0.3) +
    geom_line(aes(y = mean_cost), color = "darkgreen") +
    labs(title = "Cumulative Costs Over Time",
         y = "Total Cost") +
    theme_minimal() +
    scale_y_continuous(labels = scales::dollar_format())
  
  # Combine plots
  combined_plot <- (p1 + p2) / (p3 + p4) +
    plot_annotation(title = paste(country_name, "- Disease Dynamics Summary"),
                    theme = theme(plot.title = element_text(hjust = 0.5)))
  
  # Save plot
  ggsave(paste0(country_name, "_summary_plot.png"), 
         combined_plot, 
         width = 15, 
         height = 10)
  
  return(combined_plot)
}

# Get list of country names
country_names <- unique(Input_parameters$Country)

# Create plots for each country
for(country in country_names) {
  create_country_summary_plot(country)
}



#create country comparison violin plot monthly infections
create_country_comparison_violin <- function(country_names) {
  library(ggplot2)
  library(dplyr)
  
  # Initialize empty list to store each country's data
  all_country_data <- list()
  
  # Load and process data for each country
  for(country in country_names) {
    populations <- read.csv(paste0(country, "_populations.csv"))
    
    # Get monthly infections data, excluding first month
    country_data <- populations %>%
      filter(Month >= 2) %>%
      select(Month, Total_New_Exposures) %>%
      mutate(Country = country)  # Changed this line to use mutate
    
    all_country_data[[country]] <- country_data
  }
  
  # Combine all data
  combined_data <- do.call(rbind, all_country_data)
  
  # Create violin plot
  ggplot(combined_data, aes(x = Country, y = Total_New_Exposures)) +
    geom_violin(aes(fill = Country), alpha = 0.5) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position = "none"
    ) +
    labs(
      title = "Distribution of Monthly Infections by Country",
      x = "",
      y = "Monthly Infections"
    ) +
    scale_y_continuous(labels = scales::comma) +
    coord_cartesian(clip = "off")
}

# Use the function
country_names <- unique(Input_parameters$Country)
violin_plot <- create_country_comparison_violin(country_names)
ggsave("country_infections_violin.png", violin_plot, width = 12, height = 8)






#country exposure and infection comparison violin plots

create_country_comparison_violins <- function(country_names) {
  library(ggplot2)
  library(dplyr)
  
  # Initialize empty list to store each country's data
  all_country_data <- list()
  
  # Load and process data for each country
  for(country in country_names) {
    populations <- read.csv(paste0(country, "_populations.csv"))
    
    # Get monthly data, excluding first month
    country_data <- populations %>%
      filter(Month >= 2) %>%
      select(Month, Total_New_Exposures, Symptomatic_Cases) %>%
      mutate(Country = country)
    
    all_country_data[[country]] <- country_data
  }
  
  # Combine all data
  combined_data <- do.call(rbind, all_country_data)
  
  # Common theme elements
  plot_theme <- theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      text = element_text(size = 12)
    )
  
  # Create violin plot for exposures
  exposure_plot <- ggplot(combined_data, aes(x = Country, y = Total_New_Exposures)) +
    geom_violin(aes(fill = Country), alpha = 0.5) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +
    plot_theme +
    labs(
      title = "Distribution of Monthly Exposures by Country",
      x = "",
      y = "Monthly Exposures (log scale)"
    ) +
    scale_y_log10(labels = scales::comma) +
    coord_cartesian(clip = "off")
  
  # Create violin plot for symptomatic cases
  cases_plot <- ggplot(combined_data, aes(x = Country, y = Symptomatic_Cases)) +
    geom_violin(aes(fill = Country), alpha = 0.5) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +
    plot_theme +
    labs(
      title = "Distribution of Monthly Symptomatic Cases by Country",
      x = "",
      y = "Monthly Symptomatic Cases (log scale)"
    ) +
    scale_y_log10(labels = scales::comma) +
    coord_cartesian(clip = "off")
  
  # Save plots
  ggsave("country_exposures_violin.png", exposure_plot, width = 12, height = 8, bg = "white")
  ggsave("country_cases_violin.png", cases_plot, width = 12, height = 8, bg = "white")
  
  return(list(exposures = exposure_plot, cases = cases_plot))
}

# Use the function
country_names <- unique(Input_parameters$Country)
violin_plots <- create_country_comparison_violins(country_names)






create_country_comparison_violins <- function(country_names) {
  library(ggplot2)
  library(dplyr)
  
  # Initialize empty list to store each country's data
  all_country_data <- list()
  
  # Load and process data for each country
  for(country in country_names) {
    populations <- read.csv(paste0(country, "_populations.csv"))
    
    # Get monthly data, excluding first month
    country_data <- populations %>%
      filter(Month >= 2) %>%
      select(Month, Total_New_Exposures, Symptomatic_Cases) %>%
      mutate(Country = country)
    
    all_country_data[[country]] <- country_data
  }
  
  # Combine all data
  combined_data <- do.call(rbind, all_country_data)
  
  # Calculate median values for ordering
  country_order_exposures <- combined_data %>%
    group_by(Country) %>%
    summarise(median_exposures = median(Total_New_Exposures)) %>%
    arrange(desc(median_exposures)) %>%
    pull(Country)
  
  country_order_cases <- combined_data %>%
    group_by(Country) %>%
    summarise(median_cases = median(Symptomatic_Cases)) %>%
    arrange(desc(median_cases)) %>%
    pull(Country)
  
  # Common theme elements
  plot_theme <- theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      text = element_text(size = 12)
    )
  
  # Create violin plot for exposures
  exposure_plot <- ggplot(combined_data, aes(x = factor(Country, levels = country_order_exposures), 
                                             y = Total_New_Exposures)) +
    geom_violin(aes(fill = Country), alpha = 0.5) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +
    plot_theme +
    labs(
      title = "Distribution of Monthly Exposures by Country",
      x = "",
      y = "Monthly Exposures (log scale)"
    ) +
    scale_y_log10(labels = scales::comma) +
    coord_cartesian(clip = "off")
  
  # Create violin plot for symptomatic cases
  cases_plot <- ggplot(combined_data, aes(x = factor(Country, levels = country_order_cases), 
                                          y = Symptomatic_Cases)) +
    geom_violin(aes(fill = Country), alpha = 0.5) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +
    plot_theme +
    labs(
      title = "Distribution of Monthly Symptomatic Cases by Country",
      x = "",
      y = "Monthly Symptomatic Cases (log scale)"
    ) +
    scale_y_log10(labels = scales::comma) +
    coord_cartesian(clip = "off")
  
  # For Exposures (in the exposure violin plot function)
  exposure_summary <- combined_data %>%
    group_by(Country) %>%
    summarise(
      # Exposures
      Mean_Exposures = mean(Total_New_Exposures),
      Median_Exposures = median(Total_New_Exposures),
      SD_Exposures = sd(Total_New_Exposures),
      Q1_Exposures = quantile(Total_New_Exposures, 0.25),
      Q3_Exposures = quantile(Total_New_Exposures, 0.75),
      Min_Exposures = min(Total_New_Exposures),
      Max_Exposures = max(Total_New_Exposures)
    ) %>%
    arrange(desc(Median_Exposures))
  
  # For Cases (in the cases violin plot function)
  cases_summary <- combined_data %>%
    group_by(Country) %>%
    summarise(
      Mean_Cases = mean(Symptomatic_Cases),
      Median_Cases = median(Symptomatic_Cases),
      SD_Cases = sd(Symptomatic_Cases),
      Q1_Cases = quantile(Symptomatic_Cases, 0.25),
      Q3_Cases = quantile(Symptomatic_Cases, 0.75),
      Min_Cases = min(Symptomatic_Cases),
      Max_Cases = max(Symptomatic_Cases)
    ) %>%
    arrange(desc(Median_Cases))
  
  # Save plots
  ggsave("country_exposures_violin.png", exposure_plot, width = 12, height = 8, bg = "white")
  ggsave("country_cases_violin.png", cases_plot, width = 12, height = 8, bg = "white")
  
  #save summary statistics
  write.csv(exposure_summary, "exposure_summary_statistics.csv", row.names = FALSE)
  write.csv(cases_summary, "cases_summary_statistics.csv", row.names = FALSE)
  
  return(list(exposures = exposure_plot, cases = cases_plot))
}

# Use the function
country_names <- unique(Input_parameters$Country)
violin_plots <- create_country_comparison_violins(country_names)


#violin plots per 100k exposures/cases


library(ggplot2)
library(dplyr)


create_country_comparison_violins_per100k <- function(country_names, Input_parameters) {
  library(ggplot2)
  library(dplyr)
  
  # Initialize empty list to store each country's data
  all_country_data <- list()
  
  # Load and process data for each country
  for(country in country_names) {
    populations <- read.csv(paste0(country, "_populations.csv"))
    
    # Get population size for this country
    country_pop <- Input_parameters$initial_population_human[Input_parameters$Country == country]
    
    # Get monthly data, excluding first month and calculate per 100k rates
    country_data <- populations %>%
      filter(Month >= 2) %>%
      select(Month, Total_New_Exposures, Symptomatic_Cases) %>%
      mutate(
        Country = country,
        Exposures_per100k = (Total_New_Exposures / country_pop) * 100000,
        Cases_per100k = (Symptomatic_Cases / country_pop) * 100000
      )
    
    all_country_data[[country]] <- country_data
  }
  
  # Combine all data
  combined_data <- do.call(rbind, all_country_data)
  
  # Calculate median values for ordering
  country_order_exposures <- combined_data %>%
    group_by(Country) %>%
    summarise(median_exposures = median(Exposures_per100k)) %>%
    arrange(desc(median_exposures)) %>%
    pull(Country)
  
  country_order_cases <- combined_data %>%
    group_by(Country) %>%
    summarise(median_cases = median(Cases_per100k)) %>%
    arrange(desc(median_cases)) %>%
    pull(Country)
  
  # Common theme elements
  plot_theme <- theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      text = element_text(size = 12)
    )
  
  # Create violin plot for exposures per 100k
  exposure_plot <- ggplot(combined_data, 
                          aes(x = factor(Country, levels = country_order_exposures), 
                              y = Exposures_per100k)) +
    geom_violin(aes(fill = Country), alpha = 0.5) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +
    plot_theme +
    labs(
      title = "Distribution of Monthly Exposures per 100,000 Population by Country",
      x = "",
      y = "Monthly Exposures per 100,000 Population"
    ) +
    scale_y_continuous(labels = scales::comma) +
    coord_cartesian(clip = "off")
  
  # Create violin plot for symptomatic cases per 100k
  cases_plot <- ggplot(combined_data, 
                       aes(x = factor(Country, levels = country_order_cases), 
                           y = Cases_per100k)) +
    geom_violin(aes(fill = Country), alpha = 0.5) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +
    plot_theme +
    labs(
      title = "Distribution of Monthly Symptomatic Cases per 100,000 Population by Country",
      x = "",
      y = "Monthly Symptomatic Cases per 100,000 Population"
    ) +
    scale_y_continuous(labels = scales::comma) +
    coord_cartesian(clip = "off")
  
  # For Exposures per 100k (in the per 100k violin plot function)
  exposure_per100k_summary <- combined_data %>%
    group_by(Country) %>%
    summarise(
      Mean_Exposures_per100k = mean(Exposures_per100k),
      Median_Exposures_per100k = median(Exposures_per100k),
      SD_Exposures_per100k = sd(Exposures_per100k),
      Q1_Exposures_per100k = quantile(Exposures_per100k, 0.25),
      Q3_Exposures_per100k = quantile(Exposures_per100k, 0.75),
      Min_Exposures_per100k = min(Exposures_per100k),
      Max_Exposures_per100k = max(Exposures_per100k)
    ) %>%
    arrange(desc(Median_Exposures_per100k))
  
  # For Cases per 100k (in the per 100k violin plot function)
  cases_per100k_summary <- combined_data %>%
    group_by(Country) %>%
    summarise(
      Mean_Cases_per100k = mean(Cases_per100k),
      Median_Cases_per100k = median(Cases_per100k),
      SD_Cases_per100k = sd(Cases_per100k),
      Q1_Cases_per100k = quantile(Cases_per100k, 0.25),
      Q3_Cases_per100k = quantile(Cases_per100k, 0.75),
      Min_Cases_per100k = min(Cases_per100k),
      Max_Cases_per100k = max(Cases_per100k)
    ) %>%
    arrange(desc(Median_Cases_per100k))
  
  #
  write.csv(exposure_per100k_summary, "exposure_per100k_summary_statistics.csv", row.names = FALSE)
  write.csv(cases_per100k_summary, "cases_per100k_summary_statistics.csv", row.names = FALSE)
  
  # Save plots
  ggsave("country_exposures_per100k_violin.png", exposure_plot, width = 12, height = 8, bg = "white")
  ggsave("country_cases_per100k_violin.png", cases_plot, width = 12, height = 8, bg = "white")
  
  return(list(exposures = exposure_plot, cases = cases_plot))
}

# Use the function
country_names <- unique(Input_parameters$Country)
violin_plots_per100k <- create_country_comparison_violins_per100k(country_names, Input_parameters)




#violin plots for costs

create_cost_violin_plots <- function(country_names, Input_parameters) {
  library(ggplot2)
  library(dplyr)
  
  # Initialize empty list to store each country's data
  all_country_data <- list()
  
  # Load and process data for each country
  for(country in country_names) {
    populations <- read.csv(paste0(country, "_populations.csv"))
    
    # Get population size for this country
    country_pop <- Input_parameters$initial_population_human[Input_parameters$Country == country]
    
    # Calculate total monthly costs and per 100k rates
    country_data <- populations %>%
      filter(Month >= 2) %>%
      mutate(
        Total_Monthly_Cost = Monthly_Cost_Covert + Monthly_Cost_Ocular + Monthly_Cost_Visceral,
        Cost_per100k = (Total_Monthly_Cost / country_pop) * 100000,
        Country = country
      )
    
    all_country_data[[country]] <- country_data
  }
  
  # Combine all data
  combined_data <- do.call(rbind, all_country_data)
  
  # Calculate median values for ordering
  country_order_costs <- combined_data %>%
    group_by(Country) %>%
    summarise(median_costs = median(Total_Monthly_Cost)) %>%
    arrange(desc(median_costs)) %>%
    pull(Country)
  
  country_order_costs_per100k <- combined_data %>%
    group_by(Country) %>%
    summarise(median_costs = median(Cost_per100k)) %>%
    arrange(desc(median_costs)) %>%
    pull(Country)
  
  # Common theme elements
  plot_theme <- theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      text = element_text(size = 12)
    )
  
  # Create violin plot for total costs
  cost_plot <- ggplot(combined_data, 
                      aes(x = factor(Country, levels = country_order_costs), 
                          y = Total_Monthly_Cost)) +
    geom_violin(aes(fill = Country), alpha = 0.5) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +
    plot_theme +
    labs(
      title = "Distribution of Monthly Costs by Country",
      x = "",
      y = "Monthly Costs (USD)"
    ) +
    scale_y_log10(labels = scales::dollar_format()) +  # Using log scale for raw costs
    coord_cartesian(clip = "off")
  
  # Create violin plot for costs per 100k
  cost_per100k_plot <- ggplot(combined_data, 
                              aes(x = factor(Country, levels = country_order_costs_per100k), 
                                  y = Cost_per100k)) +
    geom_violin(aes(fill = Country), alpha = 0.5) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +
    plot_theme +
    labs(
      title = "Distribution of Monthly Costs per 100,000 Population by Country",
      x = "",
      y = "Monthly Costs per 100,000 Population (USD)"
    ) +
    scale_y_continuous(labels = scales::dollar_format()) +
    coord_cartesian(clip = "off")
  
  #costs summary statistics
  cost_summary <- combined_data %>%
    group_by(Country) %>%
    summarise(
      # Total Costs
      Mean_Cost = mean(Total_Monthly_Cost),
      Median_Cost = median(Total_Monthly_Cost),
      SD_Cost = sd(Total_Monthly_Cost),
      Q1_Cost = quantile(Total_Monthly_Cost, 0.25),
      Q3_Cost = quantile(Total_Monthly_Cost, 0.75),
      Min_Cost = min(Total_Monthly_Cost),
      Max_Cost = max(Total_Monthly_Cost),
      # Costs per 100k
      Mean_Cost_per100k = mean(Cost_per100k),
      Median_Cost_per100k = median(Cost_per100k),
      SD_Cost_per100k = sd(Cost_per100k),
      Q1_Cost_per100k = quantile(Cost_per100k, 0.25),
      Q3_Cost_per100k = quantile(Cost_per100k, 0.75),
      Min_Cost_per100k = min(Cost_per100k),
      Max_Cost_per100k = max(Cost_per100k)
    ) %>%
    arrange(desc(Median_Cost))
  
  write.csv(cost_summary, "cost_summary_statistics.csv", row.names = FALSE)
  
  # Save plots
  ggsave("country_costs_violin.png", cost_plot, width = 12, height = 8, bg = "white")
  ggsave("country_costs_per100k_violin.png", cost_per100k_plot, width = 12, height = 8, bg = "white")
  
  return(list(costs = cost_plot, costs_per100k = cost_per100k_plot))
}

# Use the function
country_names <- unique(Input_parameters$Country)
cost_violin_plots <- create_cost_violin_plots(country_names, Input_parameters)


############################################
############  DALY calculation #############
############################################


calculate_and_plot_DALYs <- function(country_names, Input_parameters) {
  library(ggplot2)
  library(dplyr)
  
  # Constants for DALY calculations
  DW_MILD <- 0.006
  DW_VISCERAL <- 0.011
  DW_OCULAR <- 0.017
  DURATION_MILD <- 0.15
  DURATION_VISCERAL <- 0.278
  
  all_country_data <- list()
  
  for(country in country_names) {
    # Read population file
    populations <- read.csv(paste0(country, "_populations.csv"))
    
    # Get population size for this country
    country_pop <- Input_parameters$initial_population_human[Input_parameters$Country == country]
    
    # Calculate yearly totals and DALYs
    yearly_data <- populations %>%
      filter(Month >= 2) %>%  # Skip first month
      mutate(Year = ceiling(Month/12)) %>%
      group_by(Year, Iteration) %>%
      summarise(
        Mild_Cases = sum(Covert_Cases),
        Visceral_Cases = sum(Visceral_Cases),
        Ocular_Cases = sum(Ocular_Cases),
        .groups = 'drop'
      ) %>%
      mutate(
        DALYs_Mild = Mild_Cases * DW_MILD * DURATION_MILD,
        DALYs_Visceral = Visceral_Cases * DW_VISCERAL * DURATION_VISCERAL,
        DALYs_Ocular = Ocular_Cases * DW_OCULAR,
        Total_DALYs = DALYs_Mild + DALYs_Visceral + DALYs_Ocular,
        DALYs_per100k = (Total_DALYs / country_pop) * 100000,
        Country = country
      )
    
    all_country_data[[country]] <- yearly_data
  }
  
  # Combine all data
  combined_data <- do.call(rbind, all_country_data)
  
  # Calculate ordering by median values
  country_order_dalys <- combined_data %>%
    group_by(Country) %>%
    summarise(median_dalys = median(Total_DALYs)) %>%
    arrange(desc(median_dalys)) %>%
    pull(Country)
  
  country_order_dalys_per100k <- combined_data %>%
    group_by(Country) %>%
    summarise(median_dalys = median(DALYs_per100k)) %>%
    arrange(desc(median_dalys)) %>%
    pull(Country)
  
  # Common theme
  plot_theme <- theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      text = element_text(size = 12)
    )
  
  # Create violin plot for total DALYs
  dalys_plot <- ggplot(combined_data, 
                       aes(x = factor(Country, levels = country_order_dalys), 
                           y = Total_DALYs)) +
    geom_violin(aes(fill = Country), alpha = 0.5) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +
    plot_theme +
    labs(
      title = "Distribution of Yearly DALYs by Country",
      x = "",
      y = "DALYs per Year"
    ) +
    scale_y_log10(labels = scales::comma) +
    coord_cartesian(clip = "off")
  
  # Create violin plot for DALYs per 100k
  dalys_per100k_plot <- ggplot(combined_data, 
                               aes(x = factor(Country, levels = country_order_dalys_per100k), 
                                   y = DALYs_per100k)) +
    geom_violin(aes(fill = Country), alpha = 0.5) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +
    plot_theme +
    labs(
      title = "Distribution of Yearly DALYs per 100,000 Population by Country",
      x = "",
      y = "DALYs per 100,000 Population"
    ) +
    scale_y_continuous(labels = scales::comma) +
    coord_cartesian(clip = "off")
  
  # Create summary statistics
  daly_summary <- combined_data %>%
    group_by(Country) %>%
    summarise(
      Mean_DALYs = mean(Total_DALYs),
      Median_DALYs = median(Total_DALYs),
      SD_DALYs = sd(Total_DALYs),
      Q1_DALYs = quantile(Total_DALYs, 0.25),
      Q3_DALYs = quantile(Total_DALYs, 0.75),
      Min_DALYs = min(Total_DALYs),
      Max_DALYs = max(Total_DALYs),
      Mean_DALYs_per100k = mean(DALYs_per100k),
      Median_DALYs_per100k = median(DALYs_per100k),
      SD_DALYs_per100k = sd(DALYs_per100k),
      Q1_DALYs_per100k = quantile(DALYs_per100k, 0.25),
      Q3_DALYs_per100k = quantile(DALYs_per100k, 0.75),
      Min_DALYs_per100k = min(DALYs_per100k),
      Max_DALYs_per100k = max(DALYs_per100k)
    ) %>%
    arrange(desc(Median_DALYs))  # Order by median DALYs
  
  # Save to CSV
  write.csv(daly_summary, "daly_summary_statistics.csv", row.names = FALSE)
  
  # Save plots
  ggsave("yearly_dalys_violin.png", dalys_plot, width = 12, height = 8, bg = "white")
  ggsave("yearly_dalys_per100k_violin.png", dalys_per100k_plot, width = 12, height = 8, bg = "white")
  
  return(list(
    dalys = dalys_plot, 
    dalys_per100k = dalys_per100k_plot,
    data = combined_data
  ))
}

# Use the function
daly_results <- calculate_and_plot_DALYs(country_names, Input_parameters)

################################################
#####Comparison plots###########################
################################################


# Read in summary statistics and calculate annual values
model_stats <- read.csv("cases_summary_statistics.csv") %>%
  mutate(
    Model_Mean = Mean_Cases * 12,
    Model_Min = Min_Cases * 12,
    Model_Max = Max_Cases * 12
  ) %>%
  select(Country, Model_Mean, Model_Min, Model_Max)

# Create comparison data by joining with input parameters
comparison_data <- data.frame(
  Country = Input_parameters$Country,
  'DALY Mean' = Input_parameters$estimated_infections,
  'DALY Min' = Input_parameters$estimated_infections_lowCI,
  'DALY Max' = Input_parameters$estimated_infections_uppCI
) %>%
  left_join(model_stats, by = "Country") %>%
  rename(
    'Model Mean' = Model_Mean,
    'Model Min' = Model_Min,
    'Model Max' = Model_Max
  )

create_comparison_plot <- function(df) {
  # First clean up the column names
  names(df) <- gsub("\\.", " ", names(df))
  
  # Reshape data to long format
  plot_data <- df %>%
    mutate(across(-Country, as.numeric)) %>%
    mutate(across(-Country, ~ifelse(. <= 0, NA, .))) %>%
    pivot_longer(
      cols = -Country,
      names_to = c("Type", "Metric"),
      names_sep = " ",
      values_to = "Value"
    ) %>%
    pivot_wider(
      names_from = Metric,
      values_from = Value
    ) %>%
    mutate(Type = factor(Type, levels = c("DALY", "Model")))
  
  # Order countries by Estimated (DALY) values
  country_order <- plot_data %>%
    filter(Type == "DALY") %>%
    arrange(desc(Mean)) %>%
    pull(Country)
  
  # Create plot
  ggplot(plot_data, aes(x = factor(Country, levels = country_order), 
                        fill = Type)) +
    geom_boxplot(
      aes(
        middle = Mean,
        lower = Min,
        upper = Max,
        ymin = Min,
        ymax = Max
      ),
      stat = "identity",
      position = position_dodge(width = 0.8),
      width = 0.6,
      alpha = 0.7
    ) +
    scale_fill_manual(values = c("DALY" = "#2171b5", "Model" = "#ef3b2c"),
                      labels = c("Input Estimate", "Model Prediction")) +
    scale_y_log10(
      labels = scales::comma,
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      minor_breaks = NULL
    ) +
    annotation_logticks(sides = "l") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      text = element_text(size = 12)
    ) +
    labs(
      title = "Comparison of Input Estimates vs Model Predictions",
      subtitle = "Note: Zero values are excluded from the log scale",
      y = "Annual Symptomatic Cases (log scale)",
      x = "",
      fill = "Data Source"
    )
}

# Create and save plot
comparison_plot <- create_comparison_plot(comparison_data)
ggsave("case_validation_comparison.png", comparison_plot, width = 12, height = 8, bg = "white")
write.csv(comparison_data, "comparison_data.csv", row.names = FALSE)

#################################################
##############Lin's correlation plots############
#################################################

library(epiR)
library(ggplot2)

# First create the ALL countries analysis and plot
analysis_data_all <- comparison_data %>%
  select(Country, DALY.Mean, `Model Mean`) %>%
  rename(Estimated = DALY.Mean,
         Predicted = `Model Mean`)

# Calculate Lin's CCC for all countries
ccc_result_all <- epi.ccc(analysis_data_all$Estimated, 
                          analysis_data_all$Predicted,
                          ci = "z-transform",
                          conf.level = 0.95)

# Create dataframe of results for all countries
ccc_summary_all <- data.frame(
  Metric = c("CCC", "95% CI Lower", "95% CI Upper", 
             "Pearson correlation (precision)", "Bias correction factor (accuracy)",
             "Scale shift", "Location shift relative to scale"),
  Value = c(ccc_result_all$rho.c[1], ccc_result_all$rho.c[2], ccc_result_all$rho.c[3],
            ccc_result_all$r[1], ccc_result_all$C.b,
            ccc_result_all$s.shift, ccc_result_all$l.shift)
)

# Save results for all countries
write.csv(ccc_summary_all, "concordance_statistics_all_countries.csv", row.names = FALSE)

# Create scatter plot for all countries
plot_all <- ggplot(analysis_data_all, aes(x = Estimated, y = Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_text(aes(label = Country), hjust = -0.2, vjust = -0.2, size = 3) +  # Added labels
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  theme_minimal() +
  labs(
    title = "Model Predictions vs Input Estimates (All Countries)",
    subtitle = paste("Lin's CCC =", round(ccc_result_all$rho.c[1], 4),
                     "\n95% CI:", round(ccc_result_all$rho.c[2], 4), "to", round(ccc_result_all$rho.c[3], 4)),
    x = "Input Estimates (log scale)",
    y = "Model Predictions (log scale)"
  ) +
  theme(
    text = element_text(size = 12),
    plot.subtitle = element_text(size = 10)
  )

# Save all countries plot
ggsave("concordance_plot_all_countries.png", plot_all, width = 10, height = 8, bg = "white")

# Now do the VALIDATION countries analysis
calibration_countries <- c("Argentina", "Brazil", "China", "Egypt", "Malaysia", 
                           "Mexico", "Nigeria", "Poland", "Romania", "Thailand", "Turkey")

analysis_data_val <- comparison_data %>%
  filter(!Country %in% calibration_countries) %>%
  select(Country, DALY.Mean, `Model Mean`) %>%
  rename(Estimated = DALY.Mean,
         Predicted = `Model Mean`)

# Calculate Lin's CCC for validation countries
ccc_result_val <- epi.ccc(analysis_data_val$Estimated, 
                          analysis_data_val$Predicted,
                          ci = "z-transform",
                          conf.level = 0.95)

# Create dataframe of results for validation countries
ccc_summary_val <- data.frame(
  Metric = c("CCC", "95% CI Lower", "95% CI Upper", 
             "Pearson correlation (precision)", "Bias correction factor (accuracy)",
             "Scale shift", "Location shift relative to scale"),
  Value = c(ccc_result_val$rho.c[1], ccc_result_val$rho.c[2], ccc_result_val$rho.c[3],
            ccc_result_val$r[1], ccc_result_val$C.b,
            ccc_result_val$s.shift, ccc_result_val$l.shift)
)

# Save results for validation countries
write.csv(ccc_summary_val, "concordance_statistics_validation_countries.csv", row.names = FALSE)

# Create scatter plot for validation countries
plot_val <- ggplot(analysis_data_val, aes(x = Estimated, y = Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_text(aes(label = Country), hjust = -0.2, vjust = -0.2, size = 3) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  theme_minimal() +
  labs(
    title = "Model Predictions vs Input Estimates (Validation Countries Only)",
    subtitle = paste("Lin's CCC =", round(ccc_result_val$rho.c[1], 4),
                     "\n95% CI:", round(ccc_result_val$rho.c[2], 4), "to", round(ccc_result_val$rho.c[3], 4)),
    x = "Input Estimates (log scale)",
    y = "Model Predictions (log scale)"
  ) +
  theme(
    text = element_text(size = 12),
    plot.subtitle = element_text(size = 10)
  )

# Save validation countries plot
ggsave("concordance_plot_validation_countries.png", plot_val, width = 10, height = 8, bg = "white")



library(colorspace)
library(epiR)
library(ggplot2)
library(ggrepel)

# Function to calculate CCC for different combinations
calculate_range_ccc <- function(data) {
  # Best case scenario (Max Model vs Max Input)
  best_case <- epi.ccc(data$DALY.Max, data$`Model Max`,
                       ci = "z-transform", conf.level = 0.95)
  
  # Worst case scenario (Min Model vs Min Input)
  worst_case <- epi.ccc(data$DALY.Min, data$`Model Min`,
                        ci = "z-transform", conf.level = 0.95)
  
  # Mean/central estimate (as before)
  mean_case <- epi.ccc(data$DALY.Mean, data$`Model Mean`,
                       ci = "z-transform", conf.level = 0.95)
  
  return(list(mean = mean_case, best = best_case, worst = worst_case))
}

create_range_plot <- function(data, title, analysis_results) {
  # Create a more diverse color palette using distinct_colors from colorspace
  n_countries <- nrow(data)
  country_colors <- colorspace::qualitative_hcl(n_countries, palette = "Dark3", 
                                                c = 80, l = 60)
  
  # Add colors to the data
  data$country_color <- country_colors
  
  ggplot(data) +
    # Add error bars for model predictions
    geom_errorbar(aes(x = DALY.Mean, ymin = `Model Min`, ymax = `Model Max`,
                      color = Country),
                  width = 0.1, alpha = 0.7) +
    # Add error bars for input estimates
    geom_errorbar(aes(xmin = DALY.Min, xmax = DALY.Max, y = `Model Mean`,
                      color = Country),
                  width = 0.1, alpha = 0.7) +
    # Add points for mean values
    geom_point(aes(x = DALY.Mean, y = `Model Mean`, color = Country)) +
    # Add line of perfect concordance
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    # Add country labels with ggrepel in black
    ggrepel::geom_text_repel(
      aes(x = DALY.Mean, y = `Model Mean`, label = Country),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      force = 2,
      segment.alpha = 0.5,
      color = "black"  # Fixed black color for labels
    ) +
    # Scale and theme settings
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(values = country_colors) +
    theme_minimal() +
    labs(
      title = title,
      subtitle = paste("Mean CCC =", round(analysis_results$mean$rho.c[1], 4),
                       "\nRange:", round(analysis_results$worst$rho.c[1], 4), 
                       "to", round(analysis_results$best$rho.c[1], 4)),
      x = "Input Estimates (log scale)",
      y = "Model Predictions (log scale)"
    ) +
    theme(
      text = element_text(size = 12),
      plot.subtitle = element_text(size = 10),
      legend.position = "none"
    )
}

# Calculate CCCs for all countries
all_results <- calculate_range_ccc(comparison_data)

# Create and save summary statistics
all_summary <- data.frame(
  Scenario = c("Mean Estimate", "Best Case", "Worst Case"),
  CCC = c(all_results$mean$rho.c[1], 
          all_results$best$rho.c[1], 
          all_results$worst$rho.c[1]),
  CI_Lower = c(all_results$mean$rho.c[2], 
               all_results$best$rho.c[2], 
               all_results$worst$rho.c[2]),
  CI_Upper = c(all_results$mean$rho.c[3], 
               all_results$best$rho.c[3], 
               all_results$worst$rho.c[3])
)
write.csv(all_summary, "concordance_statistics_all_countries_with_ranges4.csv", row.names = FALSE)

# Create plot for all countries
plot_all <- create_range_plot(comparison_data, 
                              "Model Predictions vs Input Estimates (All Countries)",
                              all_results)
ggsave("concordance_plot_all_countries_with_ranges3.png", plot_all, 
       width = 10, height = 8, bg = "white")

# Repeat for validation countries
validation_data <- comparison_data %>%
  filter(!Country %in% c("Argentina", "Brazil", "China", "Egypt", "Malaysia", 
                         "Mexico", "Nigeria", "Poland", "Romania", "Thailand", "Turkey"))

validation_results <- calculate_range_ccc(validation_data)

validation_summary <- data.frame(
  Scenario = c("Mean Estimate", "Best Case", "Worst Case"),
  CCC = c(validation_results$mean$rho.c[1], 
          validation_results$best$rho.c[1], 
          validation_results$worst$rho.c[1]),
  CI_Lower = c(validation_results$mean$rho.c[2], 
               validation_results$best$rho.c[2], 
               validation_results$worst$rho.c[2]),
  CI_Upper = c(validation_results$mean$rho.c[3], 
               validation_results$best$rho.c[3], 
               validation_results$worst$rho.c[3])
)
write.csv(validation_summary, "concordance_statistics_validation_countries_with_ranges4.csv", row.names = FALSE)

plot_validation <- create_range_plot(validation_data, 
                                     "Model Predictions vs Input Estimates (Validation Countries Only)",
                                     validation_results)
ggsave("concordance_plot_validation_countries_with_ranges3.png", plot_validation, 
       width = 10, height = 8, bg = "white")