ggplot(averaged_testing_matrix, aes(x = factor(Dose), y = CatBoost)) +
  stat_summary(aes(y = Dose, fill = "Real Dose"), 
               fun = mean, 
               geom = "bar", 
               width = 0.6, 
               alpha = 1,
               color = "steelblue4",
               linewidth = 0.8) +
  geom_jitter(aes(color = "Predicted Dose"),
              width = 0.15, 
              height = 0, 
              size = 2.5, 
              alpha = 0.8) +
  scale_fill_manual(name = NULL, 
                    values = c("Real Dose" = "honeydew3")) +
  scale_color_manual(name = NULL, 
                     values = c("Predicted Dose" = "darkorange3")) +
  labs(title = "Predicted Values by Dose Group",
       x = "Dose (Gy)",
       y = "Predicted Dose (Gy)") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10),
    legend.position = "top"
  )