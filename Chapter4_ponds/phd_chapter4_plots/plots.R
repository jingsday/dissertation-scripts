library(ggplot2)

x_range <- c(-3.5, 3.5)
y_range <- c(-1, 1)  # Similarly for y-range if needed

#Paritial plots
library(visreg)
common_theme <- theme(
  axis.title = element_text(size = 40, face = 'bold'),  # Set axis title size and bold
  axis.text = element_text(size = 36),                 # Set axis text size
  panel.background = element_rect(fill = "white"),     # Set background color to white
  panel.grid.major = element_line(color = "grey"),     # Set grid lines to grey
  panel.grid.minor = element_line(color = "lightgrey") # Optional: minor grid lines
)

df_modelling$Strat <- factor(df_modelling$Strat, levels = c(0, 1), labels = c("Absence", "Existence"))

TN_glmm_step_depth <- lmer(
  formula = log_TN ~ season_col + Cropland_500.s + Pond_dries_cat + Depth.s + sum_P.s 
  + Strat + Animals_cont.s + (1 | PondCode) + (1 | Country:Pondscape) + (1 | Country) 
  + season_col:Pond_dries_cat, data=df_modelling, na.action = "na.fail")

tab_model(TN_glmm_step_depth,transform=NULL,show.intercept = TRUE)
#interactions only season and pond dries
plots_TN <- visreg(TN_glmm_step_depth, gg = TRUE,partial=TRUE,scale='response')
plots_TN <- plots_TN[c(2, 4, 5, 7)]
plots_TN


plot_labels <- list(
  list(y = "TN\nf(Cropland 500m)", x = "Cropland (500m)"),
  list(y = "TN\nf(Depth)", x = "Depth"),
  list(y = "TN\nf(Annual precipitation)", x = "Annual precipitation"),
  list(y = "TN\nf(Livestock index)", x = "Livestock index")
)

plots_TN_minimal <- lapply(seq_along(plots_TN), function(i) {
  plots_TN[[i]] + 
    common_theme +  # Add common theme
    coord_cartesian(xlim = x_range, ylim = y_range)+# Standardize x and y axis ranges
    labs(  # Add specific labels for each plot
      y = plot_labels[[i]]$y,
      x = plot_labels[[i]]$x
    )
})

#Strat plot
vis <- visreg(TN_glmm_step_depth, "Strat", partial = TRUE, plot = FALSE)
vis_data <- vis$res  # contains partial residuals: fit, visregFit, visregRes, season_col


strat_plot_TN <- ggplot(df_modelling, aes(x = Strat, y = log_TN)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(data = vis_data, aes(x = Strat, y = visregRes), 
              width = 0.1, alpha = 0.5,size= 1) +
  theme_minimal() +
  labs(y = "TN\nf(Water Stratification)", x = "Water Stratification") 
  
strat_plot_TN <- strat_plot_TN+common_theme
strat_plot_TN


#Interactions
eff <- Effect(c("season_col", "Pond_dries_cat"), TN_glmm_step_depth, se = TRUE)

eff_df <- as.data.frame(eff)

vis <- visreg(TN_glmm_step_depth, "season_col", by = "Pond_dries_cat",partial=T, gg = TRUE, scale = "response")
vis_data <- vis$res  # contains partial residuals: fit, visregFit, visregRes, season_col


season_dryplot_TN <- ggplot(eff_df, aes(x = season_col, y = log_TN, fill = Pond_dries_cat)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(data = vis_data, aes(x = season_col, y = visregRes), 
              width = 0.1, alpha = 0.5,size= 1) +
  theme_minimal() +
  labs(y = "TN\nf(Water Stratification)", x = "Water Stratification",fill = "Hydroperiod regime") 
season_dryplot_TN


season_pond_TN <- season_pond_TN + common_theme
season_pond_TN

season_pond_TN <- ggplot(eff_df, aes(x = season_col, y = fit, fill = Pond_dries_cat)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(data = vis_data, aes(x = season_col, y = visregRes), 
              width = 0.1, alpha = 0.5,size= 1) +
  theme_minimal() +
  labs(y = "TN\nf(Water Stratification)", x = "Water Stratification",fill = "Pond category") 


season_pond_TN

list(y = "TN\nf(Hydroperiod Regime)", x = "Season")


plots_TN_minimal
TP_glmm_step_area <- lmer(log_TP ~ season_col + Cropland_500.s + Area.s + Depth.s +
                            Animals_cont.s + (1 | PondCode) + (1 | Country:Pondscape) +
                            (1 | Country) + season_col:Area.s, data = df_modelling, na.action = "na.fail")


plots_TP <- visreg(TP_glmm_step_area, gg = TRUE,partial=TRUE,scale='response')
plots_TP <- plots_TP[c(2,4,3, 5)]
plots_TP

# Define labels for each plot
plot_labels <- list(
  list(y = "TP\nf(Cropland 500m)", x = "Cropland (500m)"),
  list(y = "TP\nf(Depth)", x = "Depth"),
  list(y = "TP\nf(Area)", x = "Area"),
  list(y = "TP\nf(Livestock index)", x = "Livestock index")
)


x_range <- c(-3.5, 3.5)
y_range <- c(-2, 0.5) 

plots_TP_minimal <- lapply(seq_along(plots_TP), function(i) {
  plots_TP[[i]] + 
    common_theme +  # Add common theme
    coord_cartesian(xlim = x_range, ylim = y_range) +  # Standardize x and y axis ranges
    labs(  # Add specific labels for each plot
      y = plot_labels[[i]]$y,
      x = plot_labels[[i]]$x
    )
})




all_plots <- c( plots_TN_minimal,plots_TP_minimal)
all_plots <- c(plots_TN_minimal,list(season_pond_TN), plots_TP_minimal )

all_plots

full_combined_plot <- wrap_plots(all_plots,ncol=4 )
full_combined_plot <- full_combined_plot +
  theme(plot.margin = margin(30, 30, 30, 30)) 

# Wrap it like this to keep layout tidy
all_plots <- c(plots_TN_minimal[1:2], list(season_pond_TN), plots_TN_minimal[3:4])
full_combined_plot <- wrap_plots(all_plots, ncol = 4)

# Save the combined plot as a PNG file
out_dir <- '~/Desktop/check.png'
ggsave(out_dir, plot = full_combined_plot, width = 33, height = 22, dpi = 300)


library(patchwork)

all_plots <- c(plots_TN_minimal, plots_TP_minimal, list(season_pond_TN))
full_combined_plot <- wrap_plots(all_plots, ncol = 4)
full_combined_plot <- full_combined_plot + 
  theme(plot.margin = margin(30, 30, 30, 30))
full_combined_plot

library(effects)
library(ggplot2)

eff <- Effect(c("season_col", "Pond_dries_cat"), TN_glmm_step_depth, se = TRUE)
eff_df <- as.data.frame(eff)

# Boxplot-style plot using predicted values (not actual residuals!)
season_pond_TN <- ggplot(eff_df, aes(x = season_col, y = fit, fill = Pond_dries_cat)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.8), width = 0.2) +
  labs(
    x = "Season",
    y = "Predicted TN",
    fill = "Pond category"
  ) +
  common_theme +
  theme(
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 32),
    legend.position = "right",        # Or use c(0.8, 0.7) to place it manually
    legend.direction = "vertical"     # Make it vertical
  )

season_pond_TN
