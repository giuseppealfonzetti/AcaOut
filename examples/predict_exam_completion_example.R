# Example analysis for exam-level completion probabilities
#
# Run from the package root, e.g.:
# Rscript inst/examples/predict_exam_completion_example.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(AcaOut)
})

# Load fitted model produced in experiments/prev_fit
fit <- readRDS("experiments/prev_fit/fullBFGS.rds")

cutoff_year <- 2L
horizon <- 1L
year_length <- 365
grad_extension <- 1.5

# Exam-level completion probabilities within one year after the cutoff
exam_pred <- predict_exam_completion(
  fit,
  CUTOFF_YEAR = cutoff_year,
  HORIZON = horizon,
  YEAR_LENGTH = year_length,
  GRAD_EXTENSION = grad_extension
)

# Actual completions between cutoff and target year
cutoff_day <- cutoff_year * year_length
target_day <- as.integer(round(
  (cutoff_year + horizon - 1L) * year_length + grad_extension * year_length
))

time_long <- fit$data$timeMat |>
  as_tibble(.name_repair = "minimal") |>
  mutate(subject_id = fit$data$labs$obs) |>
  relocate(subject_id) |>
  pivot_longer(
    cols = -subject_id,
    names_to = "exam",
    values_to = "time"
  )

actual_window <- time_long |>
  semi_join(exam_pred, by = c("subject_id", "exam")) |>
  mutate(
    completed_in_window = !is.na(time) &
      time > cutoff_day &
      time <= target_day
  ) |>
  select(subject_id, exam, completed_in_window)

exam_eval <- exam_pred |>
  left_join(actual_window, by = c("subject_id", "exam")) |>
  filter(completed_in_window) |>
  group_by(exam) |>
  summarise(
    n_completed = n(),
    mean_pred_prob = mean(predicted_completion_prob),
    median_pred_prob = median(predicted_completion_prob),
    share_below_half = mean(predicted_completion_prob < 0.5),
    .groups = "drop"
  ) |>
  arrange(mean_pred_prob)

exam_eval |>
  print(n = 100)

# Outcome-level accuracy split by whether EC0101 and EC0074 were passed by cutoff
pred_one_year <- predict_k_years(
  fit,
  CUTOFF_YEAR = cutoff_year,
  HORIZON = horizon,
  YEAR_LENGTH = year_length,
  GRAD_EXTENSION = grad_extension
)

actual_outcomes <- tibble(
  subject_id = names(fit$data$outcome),
  outcome = fit$data$outcome,
  last_year = fit$data$last_year
) |>
  mutate(
    actual_within_horizon = case_when(
      outcome == 1L & last_year <= cutoff_year + horizon ~ "Graduation",
      outcome == 2L & last_year <= cutoff_year + horizon ~ "Dropout",
      outcome == 3L & last_year <= cutoff_year + horizon ~ "Transfer",
      TRUE ~ "Enrolled"
    )
  )

classification <- pred_one_year |>
  left_join(actual_outcomes, by = "subject_id") |>
  mutate(
    predicted_outcome = pmap_chr(
      list(prob_graduation, prob_dropout, prob_transfer, prob_still_enrolled),
      ~ c("Graduation", "Dropout", "Transfer", "Enrolled")[which.max(c(...))]
    ),
    correct = predicted_outcome == actual_within_horizon
  )

exam_cols <- match(c("EC0101", "EC0074"), fit$data$labs$exams)

passed_flag <- tibble(subject_id = fit$data$labs$obs) |>
  mutate(
    passed_bottlenecks = map_lgl(seq_along(subject_id), function(i) {
      in_plan <- fit$data$todoMat[i, exam_cols] == 1
      if (!any(in_plan)) {
        return(FALSE)
      }
      times <- fit$data$timeMat[i, exam_cols][in_plan]
      all(!is.na(times) & times <= cutoff_day)
    })
  )

accuracy_split <- classification |>
  left_join(passed_flag, by = "subject_id") |>
  mutate(passed_bottlenecks = replace_na(passed_bottlenecks, FALSE)) |>
  group_by(passed_bottlenecks, actual_within_horizon) |>
  summarise(
    n = n(),
    accuracy = mean(correct),
    .groups = "drop"
  )

print(accuracy_split)

#### Plots ####
year_length <- 365
max_day <- max(fit$data$timeMat, na.rm = TRUE)
grid_days <- seq(180, max_day, by = 180)

todo <- fit$data$todoMat
times <- fit$data$timeMat
exam_labels <- fit$data$labs$exams

exam_info <- tibble(exam = exam_labels) |>
  mutate(
    median_day = map_dbl(seq_along(exam_labels), function(j) {
      vals <- times[, j][todo[, j] == 1]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) {
        return(NA_real_)
      }
      median(vals)
    }),
    median_year = if_else(
      is.na(median_day),
      NA_real_,
      ceiling(median_day / year_length)
    ),
    passing_curve = map(seq_along(exam_labels), function(j) {
      vals <- times[, j][todo[, j] == 1]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) {
        return(tibble(day = numeric(0), cum_rate = numeric(0)))
      }
      tibble(
        day = grid_days,
        cum_rate = map_dbl(grid_days, ~ mean(vals <= .x))
      )
    })
  )

exam_times_long <- exam_info |>
  filter(!is.na(median_year)) |>
  select(exam, median_year, passing_curve) |>
  unnest(passing_curve) |>
  mutate(
    years_since_enrolment = day / year_length,
    year_bucket = factor(median_year)
  ) |>
  mutate(
    highlight = exam %in% c("EC0101", "EC0074")
  )

plot_obj <- ggplot(
  exam_times_long,
  aes(
    x = years_since_enrolment,
    y = cum_rate,
    group = exam,
    colour = highlight
  )
) +
  geom_line() +
  facet_wrap(~year_bucket) +
  labs(
    x = "Years since enrolment",
    y = "Cumulative passing rate",
    colour = "Median completion year",
    title = "Exam-level completion curves grouped by completion year"
  ) +
  theme_minimal()
plot_obj


#######

fit <- fullBFGS

year_length <- 365
max_day <- max(fit$data$timeMat, na.rm = TRUE)
grid_days <- seq(180, max_day, by = 180)

todo <- fit$data$todoMat
times <- fit$data$timeMat
exam_labels <- fit$data$labs$exams

bottleneck_exams <- c("EC0101", "EC0074")

exam_info <- tibble(exam = exam_labels) |>
  mutate(
    median_day = map_dbl(seq_along(exam_labels), function(j) {
      vals <- times[, j][todo[, j] == 1]
      if (length(vals) == 0) {
        return(NA_real_)
      }
      median(vals, na.rm = TRUE)
    }),
    median_year = if_else(
      is.na(median_day),
      NA_real_,
      ceiling(median_day / year_length)
    ),
    passing_curve = map(seq_along(exam_labels), function(j) {
      vals <- times[, j][todo[, j] == 1]
      if (length(vals) == 0) {
        return(tibble(day = numeric(0), cum_rate = numeric(0)))
      }
      vals <- replace_na(vals, Inf)
      tibble(
        day = grid_days,
        cum_rate = map_dbl(grid_days, ~ mean(vals <= .x, na.rm = TRUE))
      )
    }),
    is_bottleneck = exam %in% bottleneck_exams
  )

exam_times_long <- exam_info |>
  filter(!is.na(median_year)) |>
  select(exam, median_year, passing_curve, is_bottleneck) |>
  unnest(passing_curve) |>
  mutate(
    years_since_enrolment = day / year_length,
    year_bucket = factor(median_year)
  )

plot_obj <- ggplot(
  exam_times_long,
  aes(
    x = years_since_enrolment,
    y = cum_rate,
    group = exam,
    colour = year_bucket
  )
) +
  geom_line(alpha = 0.3) +
  geom_line(
    data = exam_times_long |> filter(is_bottleneck),
    aes(group = exam),
    colour = "#D55E00",
    linewidth = 1
  ) +
  geom_text(
    data = exam_times_long |> filter(is_bottleneck, day == max(day)),
    aes(label = exam),
    colour = "#D55E00",
    hjust = -0.1,
    size = 3
  ) +
  facet_wrap(~year_bucket) +
  labs(
    x = "Years since enrolment",
    y = "Cumulative passing rate",
    colour = "Median completion year",
    title = "Exam-level completion curves grouped by completion year",
    subtitle = "Bottleneck exams (EC0101, EC0074) highlighted"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5)) +
  theme(legend.position = "none")

limit_year <- max(exam_times_long$years_since_enrolment)
plot_obj <- plot_obj + xlim(0, limit_year + 1)
plot_obj
