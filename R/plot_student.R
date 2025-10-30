plot_student <- function(MOD, SUBJECT_ID) {
  sub_id_out <- fullBFGS$data$outcome[sub_id]
  sub_id_out <- switch(
    sub_id_out,
    "1" = "Graduation",
    "2" = "Dropout",
    "3" = "Transfer",
    "0" = "Enrolled"
  )
  sub_id_year <- fullBFGS$data$last_year[sub_id]
  sub_id_ref <- tibble(year = sub_id_year, outcome = sub_id_out)
  gg_pred <- pred_grid |>
    unnest(pred) |>
    filter(subject_id == sub_id) |>
    filter(start < sub_id_year) |>
    pivot_longer(
      cols = starts_with("prob"),
      values_to = "prob",
      names_to = "outcome"
    ) |>
    mutate(
      outcome = map_chr(
        outcome,
        ~ case_when(
          str_detect(.x, "graduation") ~ "Graduation",
          str_detect(.x, "dropout") ~ "Dropout",
          str_detect(.x, "transfer") ~ "Transfer",
          str_detect(.x, "still") ~ "Enrolled",
          TRUE ~ .x
        )
      ),
      cutoff_year_lab = paste0("Predictions at the end\nof year ", cutoff_year)
    ) |>
    ggplot(aes(x = target_year, y = prob, col = outcome)) +
    facet_grid(cutoff_year_lab ~ .) +
    geom_rect(
      aes(xmin = 0, xmax = cutoff_year, ymin = 0, ymax = 1),
      fill = "#e7e7e7",
      col = "#e7e7e7"
    ) +
    geom_vline(
      data = sub_id_ref,
      aes(xintercept = year, col = outcome),
      linetype = "dashed",
      show.legend = FALSE
    ) +
    geom_point() +
    geom_step(linewidth = 0.8, alpha = .8) +
    geom_text(
      aes(
        x = .92,
        y = .93,
        label = paste0("Exams done: ", n_exams_done, "/", n_exams_todo)
      ),
      col = "#383838"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      panel.grid.major.x = element_blank()
    ) +
    labs(x = "Academic Year", col = "", y = "Cumulative probability") +
    scale_color_manual(values = my_pal) +
    scale_fill_manual(values = my_pal)

  map_traj <- map_data |>
    filter(year_orig <= event_year) |>
    mutate(highlight = subject_id == sub_id) |>
    pivot_longer(
      cols = c(ability, speed),
      names_to = "latent",
      values_to = "val"
    ) |>
    mutate(
      latent = recode_factor(latent, "ability" = "Ability", "speed" = "Speed")
    )
  gg_post_lat <- map_traj |>
    ggplot(aes(x = year, y = val, group = subject_id)) +
    facet_grid(latent ~ ., scales = "free") +
    theme_bw() +
    geom_line(aes(alpha = if_else(highlight, 1, 0.01))) +
    geom_point(
      data = map_traj |>
        filter(event_year == year_orig),
      aes(x = year, y = val),
      alpha = .1
    ) +
    geom_point(
      data = map_traj |> filter(highlight, event_year == year_orig),
      size = 3,
      aes(x = year, col = outcome)
    ) +
    geom_line(
      data = map_traj |>
        filter(highlight),
      size = 1,
      aes(col = outcome)
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
    ) +
    labs(x = "", y = "") +
    # scale_color_okabeito() +
    scale_color_manual(values = my_pal)
}
