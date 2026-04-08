## copy of combined_analysis_all_colonies.R with the notebook plots/models added
library(dplyr)
library(ggplot2)
library(Biostrings)
library(tidyverse)
library(data.table)
library(patchwork)
library(broom)
library(scales)

project_root = "/Users/u0174312/Documents/New project"
output_dir = file.path(project_root, "results", "06_genome_wide_design_rules_from_original_script")
plots_dir = file.path(output_dir, "plots")
tables_dir = file.path(output_dir, "tables")
inputs_dir = file.path(output_dir, "inputs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)

theme_design_rules = function(aspect_ratio = 1, legend_position = "top") {
  theme_classic(base_size = 11) +
    theme(
      legend.position = legend_position,
      axis.ticks.x = element_blank(),
      aspect.ratio = aspect_ratio
    )
}

gc_fraction = function(sequence) {
  sequence = toupper(as.character(sequence))
  ifelse(
    nchar(sequence) > 0,
    (stringr::str_count(sequence, "G") + stringr::str_count(sequence, "C")) / nchar(sequence),
    NA_real_
  )
}

gc_fraction_bin = function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x < 0.30 ~ "<30%",
    x < 0.40 ~ "30-39%",
    x < 0.50 ~ "40-49%",
    x < 0.60 ~ "50-59%",
    TRUE ~ ">=60%"
  )
}

binom_wilson = function(successes, trials, conf_level = 0.95) {
  z_value = qnorm(1 - (1 - conf_level) / 2)
  rate = ifelse(trials > 0, successes / trials, NA_real_)
  denominator = 1 + (z_value^2 / trials)
  center = (rate + (z_value^2 / (2 * trials))) / denominator
  half_width = (
    z_value *
      sqrt((rate * (1 - rate) + (z_value^2 / (4 * trials))) / trials)
  ) / denominator

  tibble(
    rate = rate,
    lower = pmax(0, center - half_width),
    upper = pmin(1, center + half_width)
  )
}

summarize_edit_rates = function(data, group_cols) {
  summary_table = data %>%
    filter(assayed_design) %>%
    group_by(across(all_of(group_cols))) %>%
    summarize(
      total_guides = n(),
      successful_guides = sum(edited_design),
      total_assayed_colonies = sum(assayed_colonies),
      edited_colonies = sum(edited_colonies),
      .groups = "drop"
    )

  design_ci = binom_wilson(summary_table$successful_guides, summary_table$total_guides)
  colony_ci = binom_wilson(summary_table$edited_colonies, summary_table$total_assayed_colonies)

  bind_cols(
    summary_table,
    design_ci %>% dplyr::rename(
      design_success_rate = rate,
      design_success_lower = lower,
      design_success_upper = upper
    ),
    colony_ci %>% dplyr::rename(
      colony_edit_rate = rate,
      colony_edit_lower = lower,
      colony_edit_upper = upper
    )
  )
}

annotate_nearest_score_sv = function(data, chr_col, pos_col, score_sv_sites) {
  split_scores = split(score_sv_sites, score_sv_sites$chr)

  nearest_hits = map2(
    data[[chr_col]],
    data[[pos_col]],
    function(chr_value, pos_value) {
      chr_scores = split_scores[[as.character(chr_value)]]

      if (is.null(chr_scores) || is.na(pos_value) || nrow(chr_scores) == 0) {
        return(tibble(
          nearest_score_sv = NA_real_,
          nearest_score_sv_pos = NA_integer_,
          nearest_score_sv_guide = NA_character_,
          nearest_score_sv_site_id = NA_character_
        ))
      }

      nearest_index = which.min(abs(chr_scores$pos - as.integer(pos_value)))

      tibble(
        nearest_score_sv = as.numeric(chr_scores$score_sv[[nearest_index]]),
        nearest_score_sv_pos = as.integer(chr_scores$pos[[nearest_index]]),
        nearest_score_sv_guide = as.character(chr_scores$guide[[nearest_index]]),
        nearest_score_sv_site_id = as.character(chr_scores$id[[nearest_index]])
      )
    }
  ) %>%
    bind_rows()

  bind_cols(data, nearest_hits)
}

plot_rate_bars = function(summary_table, x_col, x_label, title_text) {
  ggplot(summary_table, aes(x = .data[[x_col]], y = design_success_rate, fill = .data[[x_col]])) +
    geom_col(width = 0.72, color = "black", linewidth = 0.2, show.legend = FALSE) +
    geom_errorbar(aes(ymin = design_success_lower, ymax = design_success_upper), width = 0.15) +
    geom_text(aes(label = paste0("n=", total_guides)), vjust = -0.35, size = 3) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05), expand = c(0, 0)) +
    labs(x = x_label, y = "Design success rate", title = title_text) +
    theme_design_rules(aspect_ratio = 1, legend_position = "none")
}

plot_model_odds_ratios = function(model_table, title_text) {
  plot_data = model_table %>%
    filter(term != "(Intercept)") %>%
    filter(!is.na(odds_ratio), !is.na(odds_ratio_low), !is.na(odds_ratio_high)) %>%
    mutate(term_label = factor(term_label, levels = rev(term_label)))

  ggplot(plot_data, aes(
    x = odds_ratio,
    y = term_label,
    xmin = odds_ratio_low,
    xmax = odds_ratio_high,
    color = feature_group
  )) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(height = 0.18, linewidth = 0.5) +
    geom_point(size = 2.4) +
    scale_x_log10() +
    scale_color_manual(values = c(
      "PAM" = "#B56576",
      "PAM upstream base" = "#E56B6F",
      "Edit position" = "#6D597A",
      "Guide score" = "#355070",
      "GC content" = "#2A9D8F",
      "DR folding" = "#7A8C5F",
      "Homopolymers" = "#6C9A8B",
      "Variant size" = "#BC6C25",
      "Round" = "#6C757D",
      "Other" = "#666666"
    )) +
    labs(x = "Odds ratio (log scale)", y = NULL, color = NULL, title = title_text) +
    theme_design_rules(aspect_ratio = NULL, legend_position = "top")
}

format_unified_model_term = function(term) {
  formatted_term = term

  replacement_map = c(
    "pamTTTC" = "TTTC vs TTTA",
    "pamTTTG" = "TTTG vs TTTA",
    "pam_dist_bin7-11" = "7-11 nt vs 1-6 nt",
    "pam_dist_bin12-17" = "12-17 nt vs 1-6 nt",
    "pam_dist_bin18-23" = "18-23 nt vs 1-6 nt",
    "pam_dist_binPAM" = "PAM-overlap vs 1-6 nt",
    "pam_dist_bin>23" = ">23 nt vs 1-6 nt",
    "pam_upstream_ntA" = "PAM 5' A vs C",
    "pam_upstream_ntG" = "PAM 5' G vs C",
    "pam_upstream_ntT" = "PAM 5' T vs C",
    "deepcpf1_score" = "DeepCpf1 score",
    "guide_gc_20" = "Guide GC (first 20 nt)",
    "donor_wt_gc" = "WT donor GC",
    "dr_spacer_pair_anyTRUE" = "Any DR-spacer pairing",
    "roundround_2" = "Round 2 vs Round 1"
  )

  for (term_name in names(replacement_map)) {
    formatted_term = stringr::str_replace_all(formatted_term, fixed(term_name), replacement_map[[term_name]])
  }

  formatted_term = stringr::str_replace_all(formatted_term, fixed(":"), " x ")
  formatted_term
}

tidy_unified_model = function(model) {
  tidy(model, conf.int = TRUE) %>%
    mutate(
      odds_ratio = exp(estimate),
      odds_ratio_low = exp(conf.low),
      odds_ratio_high = exp(conf.high),
      feature_group = case_when(
        str_detect(term, ":") ~ "Interactions",
        str_starts(term, "pam_upstream_nt") ~ "PAM upstream base",
        str_starts(term, "pam_dist_bin") ~ "Edit position",
        str_starts(term, "pam") ~ "PAM",
        term == "deepcpf1_score" ~ "Guide score",
        term %in% c("guide_gc_20", "donor_wt_gc") ~ "GC content",
        term == "dr_spacer_pair_anyTRUE" ~ "DR folding",
        str_starts(term, "round") ~ "Round",
        TRUE ~ "Other"
      ),
      term_label = vapply(term, format_unified_model_term, character(1))
    ) %>%
    filter(term != "(Intercept)")
}

save_plot = function(plot_object, filename, width, height) {
  ggsave(filename = filename, plot = plot_object, width = width, height = height, dpi = 300)
}

summarize_rule_comparison = function(preferred_data, other_data, priority, rule, preferred_group, other_group, note = "") {
  tibble(
    priority = priority,
    rule = rule,
    preferred_group = preferred_group,
    other_group = other_group,
    preferred_designs = nrow(preferred_data),
    other_designs = nrow(other_data),
    preferred_design_success_rate = mean(preferred_data$edited_design),
    other_design_success_rate = mean(other_data$edited_design),
    preferred_colony_edit_rate = sum(preferred_data$edited_colonies) / sum(preferred_data$assayed_colonies),
    other_colony_edit_rate = sum(other_data$edited_colonies) / sum(other_data$assayed_colonies),
    note = note
  )
}

summarize_design_subset = function(data, label) {
  successes = sum(data$edited_design)
  trials = nrow(data)
  ci = binom_wilson(successes, trials)

  tibble(
    filter_label = label,
    total_guides = trials,
    successful_guides = successes,
    design_success_rate = ci$rate,
    design_success_lower = ci$lower,
    design_success_upper = ci$upper
  )
}

plot_optimal_filter_effects = function(summary_table) {
  plot_data = summary_table %>%
    mutate(
      filter_label = factor(
        filter_label,
        levels = rev(c(
          "All assayed designs",
          "DeepCpf1 top quartile",
          "PAM distance not PAM/18-23",
          "PAM not TTTG",
          "5' PAM base not C",
          "All 4 optimal filters"
        ))
      ),
      filter_group = case_when(
        filter_label == "All assayed designs" ~ "Baseline",
        filter_label == "All 4 optimal filters" ~ "Combined filter",
        TRUE ~ "Single filter"
      ),
      label = sprintf("%s   n=%d", percent(design_success_rate, accuracy = 1), total_guides)
    )

  ggplot(plot_data, aes(x = design_success_rate, y = filter_label, fill = filter_group)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.2, show.legend = TRUE) +
    geom_errorbarh(aes(xmin = design_success_lower, xmax = design_success_upper), height = 0.16) +
    geom_text(aes(label = label), hjust = -0.05, size = 3.2) +
    scale_fill_manual(values = c(
      "Baseline" = "#D9D9D9",
      "Single filter" = "#6C9A8B",
      "Combined filter" = "#355070"
    )) +
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1.12),
      expand = c(0, 0)
    ) +
    labs(
      x = "Designs edited",
      y = NULL,
      fill = NULL,
      title = "Optimal cassette filters enrich successful designs"
    ) +
    theme_design_rules(aspect_ratio = NULL, legend_position = "top")
}

plot_three_way_filter_effects = function(summary_table) {
  combo_order = summary_table %>%
    filter(filter_group == "3-way combination") %>%
    arrange(design_success_rate) %>%
    pull(filter_label)

  plot_data = summary_table %>%
    mutate(
      filter_label = factor(
        filter_label,
        levels = c("All assayed designs", combo_order, "All 4 optimal filters")
      ),
      label = sprintf("%s   n=%d", percent(design_success_rate, accuracy = 1), total_guides)
    )

  ggplot(plot_data, aes(x = design_success_rate, y = filter_label, fill = filter_group)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.2) +
    geom_errorbarh(aes(xmin = design_success_lower, xmax = design_success_upper), height = 0.16) +
    geom_text(aes(label = label), hjust = -0.05, size = 3.1) +
    scale_fill_manual(values = c(
      "Baseline" = "#D9D9D9",
      "3-way combination" = "#6C9A8B",
      "All 4 filters" = "#355070"
    )) +
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1.12),
      expand = c(0, 0)
    ) +
    labs(
      x = "Designs edited",
      y = NULL,
      fill = NULL,
      title = "Three-way optimal-filter combinations"
    ) +
    theme_design_rules(aspect_ratio = NULL, legend_position = "top")
}

plot_rule_train_test_effects = function(summary_table) {
  plot_data = summary_table %>%
    mutate(
      split = factor(split, levels = c("Train", "Test")),
      subset_label = factor(subset_label, levels = rev(c("All assayed designs", "3-rule set"))),
      fill_group = ifelse(subset_label == "All assayed designs", "Baseline", "3-rule set"),
      label = sprintf("%s   n=%d", percent(design_success_rate, accuracy = 1), total_guides)
    )

  ggplot(plot_data, aes(x = design_success_rate, y = subset_label, fill = fill_group)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.2) +
    geom_errorbarh(aes(xmin = design_success_lower, xmax = design_success_upper), height = 0.16) +
    geom_text(aes(label = label), hjust = -0.05, size = 3.2) +
    facet_wrap(~ split, nrow = 1) +
    scale_fill_manual(values = c(
      "Baseline" = "#D9D9D9",
      "3-rule set" = "#355070"
    )) +
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1.12),
      expand = c(0, 0)
    ) +
    labs(
      x = "Designs edited",
      y = NULL,
      fill = NULL,
      title = "Three-rule cassette filter retains enrichment in held-out designs",
      subtitle = "Rule = PAM distance not PAM/18-23 + PAM not TTTG + 5' PAM base not C"
    ) +
    theme_design_rules(aspect_ratio = NULL, legend_position = "top")
}

plot_selected_filter_comparison = function(summary_table) {
  plot_data = summary_table %>%
    mutate(
      filter_label = factor(
        filter_label,
        levels = rev(c(
          "All assayed designs",
          "PAM distance not PAM/18-23 + PAM not TTTG + 5' PAM base not C",
          "DeepCpf1 top quartile + PAM distance not PAM/18-23 + PAM not TTTG",
          "All 4 optimal filters"
        ))
      ),
      label = sprintf("%s   n=%d", percent(design_success_rate, accuracy = 1), total_guides)
    )

  ggplot(plot_data, aes(x = design_success_rate, y = filter_label)) +
    geom_col(width = 0.72, fill = "#E69F00", color = "black", linewidth = 0.2) +
    geom_errorbarh(aes(xmin = design_success_lower, xmax = design_success_upper), height = 0.16) +
    geom_text(aes(label = label), hjust = -0.05, size = 3.2) +
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1.12),
      expand = c(0, 0)
    ) +
    labs(
      x = "Designs edited",
      y = NULL,
      title = "Editing rate across selected cassette-filter sets"
    ) +
    theme_design_rules(aspect_ratio = NULL, legend_position = "none")
}

summarize_stacked_design_outcomes = function(data, parameter_col, parameter_pretty, parameter_group) {
  filtered_data = data %>%
    filter(assayed_design, !is.na(.data[[parameter_col]]))

  if (nrow(filtered_data) == 0) {
    return(tibble())
  }

  parameter_values = filtered_data[[parameter_col]]
  if (is.factor(parameter_values)) {
    value_levels = levels(droplevels(parameter_values))
  } else {
    value_levels = sort(unique(as.character(parameter_values)))
  }

  filtered_data %>%
    transmute(
      parameter = parameter_col,
      parameter_pretty = parameter_pretty,
      parameter_group = parameter_group,
      parameter_value = as.character(.data[[parameter_col]]),
      outcome_state = outcome_state
    ) %>%
    count(parameter, parameter_pretty, parameter_group, parameter_value, outcome_state, name = "design_count") %>%
    complete(
      parameter,
      parameter_pretty,
      parameter_group,
      parameter_value = value_levels,
      outcome_state = c("Intended edit only", "Other variants observed", "Unedited"),
      fill = list(design_count = 0L)
    ) %>%
    group_by(parameter, parameter_pretty, parameter_group, parameter_value) %>%
    mutate(
      total_designs = sum(design_count),
      fraction = ifelse(total_designs > 0, design_count / total_designs, 0)
    ) %>%
    ungroup() %>%
    mutate(
      value_order = match(parameter_value, value_levels)
    )
}

plot_stacked_design_outcomes = function(summary_table, title_text, subtitle_text = NULL, ncol = 2) {
  plot_data = summary_table %>%
    mutate(
      parameter_pretty = factor(parameter_pretty, levels = unique(parameter_pretty)),
      panel_value = paste(parameter_pretty, parameter_value, sep = "___")
    )

  panel_levels = plot_data %>%
    distinct(parameter_pretty, parameter_group, parameter_value, value_order, panel_value) %>%
    arrange(parameter_pretty, value_order) %>%
    pull(panel_value)

  plot_data = plot_data %>%
    mutate(
      panel_value = factor(panel_value, levels = panel_levels),
      outcome_state = factor(
        outcome_state,
        levels = c("Unedited", "Other variants observed", "Intended edit only")
      )
    )

  total_label_data = plot_data %>%
    distinct(parameter_pretty, panel_value, total_designs)

  other_variant_label_data = plot_data %>%
    filter(outcome_state == "Other variants observed") %>%
    distinct(parameter_pretty, panel_value, other_variant_designs = design_count)

  ggplot(plot_data, aes(x = panel_value, y = fraction, fill = outcome_state)) +
    geom_col(width = 0.82, color = "white", linewidth = 0.2) +
    geom_text(
      data = total_label_data,
      aes(x = panel_value, y = 1.05, label = paste0("n=", total_designs)),
      inherit.aes = FALSE,
      size = 2.5
    ) +
    geom_text(
      data = other_variant_label_data,
      aes(x = panel_value, y = 1.015, label = paste0("other=", other_variant_designs)),
      inherit.aes = FALSE,
      size = 2.3,
      color = "#D55E00"
    ) +
    facet_wrap(~ parameter_pretty, scales = "free_x", ncol = ncol) +
    scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1.10),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = c(
        "Intended edit only" = "#18A3A6",
        "Other variants observed" = "#D55E00",
        "Unedited" = "#D9D9D9"
      )
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL,
      y = "Fraction of assayed designs",
      fill = NULL,
      title = title_text,
      subtitle = subtitle_text
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "top",
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = grid::unit(1.1, "lines"),
      plot.margin = margin(10, 10, 18, 10)
    )
}

plot_chromosome_editing_circos = function(summary_table, chromosome_layout, bin_width_bp) {
  plot_data = summary_table %>%
    mutate(
      edit_track_ymin = 0.70,
      edit_track_ymax = 0.98,
      edit_track_padding = 0.012,
      count_track_ymin = 0.46,
      count_track_ymax = 0.66,
      count_track_padding = 0.012,
      edit_bar_ymin = edit_track_ymin + edit_track_padding,
      edit_bar_ymax = edit_bar_ymin + (edit_track_ymax - edit_track_ymin - 2 * edit_track_padding) * pmin(pmax(design_success_rate, 0), 1),
      count_bar_ymin = count_track_ymin + count_track_padding,
      count_bar_ymax = count_bar_ymin + (count_track_ymax - count_track_ymin - 2 * count_track_padding) * pmin(pmax(assayed_count_scaled, 0), 1)
    )

  total_span = max(chromosome_layout$genome_end) + first(chromosome_layout$chr_gap)

  ggplot() +
    geom_rect(
      data = chromosome_layout,
      aes(xmin = genome_start, xmax = genome_end, ymin = 0.70, ymax = 0.98),
      fill = "white",
      color = "black",
      linewidth = 0.35
    ) +
    geom_rect(
      data = chromosome_layout,
      aes(xmin = genome_start, xmax = genome_end, ymin = 0.46, ymax = 0.66),
      fill = "white",
      color = "black",
      linewidth = 0.35
    ) +
    geom_rect(
      data = plot_data %>% filter(assayed_designs > 0),
      aes(xmin = bin_plot_start, xmax = bin_plot_end, ymin = count_bar_ymin, ymax = count_bar_ymax),
      fill = "#CFCFCF",
      color = NA
    ) +
    geom_rect(
      data = plot_data %>% filter(assayed_designs > 0),
      aes(xmin = bin_plot_start, xmax = bin_plot_end, ymin = edit_bar_ymin, ymax = edit_bar_ymax),
      fill = "#18A3A6",
      color = NA
    ) +
    geom_text(
      data = chromosome_layout,
      aes(
        x = genome_mid,
        y = 1.08,
        label = chr_label,
        angle = chromosome_label_angle,
        hjust = chromosome_label_hjust
      ),
      size = 3.1,
      fontface = "bold"
    ) +
    scale_x_continuous(limits = c(0, total_span), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.14), expand = c(0, 0)) +
    coord_polar(theta = "x", clip = "off") +
    labs(
      title = "Genome-wide editing across chromosome regions",
      subtitle = paste0(
        "Outer teal bars: editing rate per ", comma(bin_width_bp / 1000), " kb window. ",
        "Inner gray bars: assayed-design count per window."
      ),
      caption = paste0(
        "Editing rate is calculated on assayed designs only. Empty windows had no assayed designs. ",
        "Assayed-design count bars are scaled to the maximum occupied window (n = ", max(plot_data$assayed_designs), ")."
      )
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5, color = "grey35"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(18, 40, 28, 40)
    )
}

## original inputs
picked_colonies = read.csv("/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/REDI/picked_colonies_20250613_annotated.csv")
picked_colonies_round_2 = read.csv("/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/REDI/second_REDI_20250917/scritps_20251216/results/colonies_with_BC1_20260212.csv")
designs = read.delim("/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/Selected_812/designs_812_20241106.txt", sep = "\t", header = TRUE)
design_annotations_extra = read.delim(file.path(project_root, "scripts", "designs_812_sv_my_variants_annotated.tsv"), sep = "\t", header = TRUE)
variant_annotations = read.delim(file.path(project_root, "results", "variants", "vars_20241106.annotated.tsv"), sep = "\t", header = TRUE)

picked_colonies_simplified = picked_colonies %>%
  select(design, BC1, colony_ID) %>%
  mutate(round = "round_1")

picked_colonies_round_2_simplified = picked_colonies_round_2 %>%
  select(design, BC1, colony_id) %>%
  dplyr::rename(colony_ID = colony_id) %>%
  mutate(round = "round_2")

picked_colonies = bind_rows(picked_colonies_simplified, picked_colonies_round_2_simplified)

designs = designs %>% filter(uid %in% picked_colonies$design)

designs_simplified = designs %>%
  select(uid, guide_id, guide_seq, guide_nopam, donor_alt, donor_wt, ref_len, alt_len, vardiff, pam_dist_var_start, pam_dist_var_end, gc_content, pam, max_A, max_T, num_T, CHR, POS, REF, ALT)

picked_colonies_simplified = picked_colonies %>%
  select(design, BC1, colony_ID, round)

colnames(picked_colonies_simplified) = c("uid", "BC1", "colony_ID", "round")

picked_colonies_annotated = left_join(picked_colonies_simplified, designs_simplified, by = "uid")

picked_colonies_annotated = picked_colonies_annotated %>%
  mutate(BC_rc = as.character(reverseComplement(DNAStringSet(BC1))))

## target coverage
coverage_round_1 = read.table("/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/Tn5/Analysis/stats/target_basecov_20260222.tsv", header = TRUE)
coverage_round_1$colony_ID = sub(".*VB_(.*?)_S.*", "\\1", coverage_round_1$sample)
coverage_round_1$round = "round_1"

coverage_round_2 = read.table("/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/Tn5/second_round_20260130/read_stats/target_coverage_20260213.tsv", header = TRUE)
coverage_round_2$colony_ID = sub(".*VB_(.*?)_S.*", "\\1", coverage_round_2$sample)
coverage_round_2$round = "round_2"

target_pos_cov = bind_rows(coverage_round_1, coverage_round_2) %>%
  select(round, colony_ID, coverage, sample)

mean_cov = read.csv("/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/Tn5/second_round_20260130/read_stats/mean_cov_20260223.csv")
mean_cov_second = read.csv("/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/Tn5/second_round_20260130/read_stats/mean_cov_20260223_second.csv")
mean_cov = rbind(mean_cov, mean_cov_second) %>%
  mutate(sample = sub("\\.bincov\\.txt$", "", file)) %>%
  select(sample, mean_cov)

picked_colonies_annotated = picked_colonies_annotated %>%
  left_join(target_pos_cov, by = c("round", "colony_ID")) %>%
  dplyr::rename(target_basecov = coverage) %>%
  left_join(mean_cov, by = "sample")

designs_with_target_cov = unique((picked_colonies_annotated %>% filter(target_basecov >= 2) %>% select(uid))$uid)

## normalized target variants
dummy_info = "AB=0;ABP=0;AC=1;AF=1;AN=1;AO=8;CIGAR=1X;DP=9;DPB=9;DPRA=0;EPP=4.09604;EPPR=5.18177;GTI=0;LEN=1;MEANALT=1;MQM=43.75;MQMR=45;NS=1;NUMALT=1;ODDS=54.7201;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=328;QR=30;RO=1;RPL=1;RPP=12.7819;RPPR=5.18177;RPR=7;RUN=1;SAF=4;SAP=3.0103;SAR=4;SRF=0;SRP=5.18177;SRR=1;TYPE=snp"
dummy_format = "GT:DP:AD:RO:QR:AO:QA:GL"
dummy_unknown = "1:9:1,8:1:30:8:328:-25.5899,0"

designs_indels = designs %>%
  filter(vardiff != 0) %>%
  mutate(QUAL = 100, FILTER = ".", ID = uid, INFO = dummy_info, FORMAT = dummy_format, unknown = dummy_unknown) %>%
  select(CHR, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, unknown)

colnames(designs_indels) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "unknown")

write.table(
  designs_indels,
  file.path(inputs_dir, "vcf_indels_20260223.vcf"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

normalized_indel_loci = read.table(
  "/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/Tn5/second_round_20260130/annotations_ref/normalized_indel_designs_20260223.vcf",
  sep = "\t",
  header = FALSE,
  comment.char = "#"
) %>%
  select(V1, V2, V3, V4, V5)

colnames(normalized_indel_loci) = c("CHR", "POS", "uid", "REF", "ALT")

designs_indels = picked_colonies_annotated %>%
  filter(vardiff != 0) %>%
  select(-CHR, -POS, -REF, -ALT)

designs_indels = left_join(designs_indels, normalized_indel_loci, by = "uid") %>%
  mutate(sample_vars = paste(sample, CHR, POS, REF, ALT, sep = "_"))

designs_snps = picked_colonies_annotated %>%
  filter(vardiff == 0) %>%
  mutate(sample_vars = paste(sample, CHR, POS, REF, ALT, sep = "_"))

designs_all = rbind(designs_snps, designs_indels)

normalized_vcf_parsed = read.csv("/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/Tn5/Analysis/results/variants_norm_WGS_20260222.csv") %>%
  select(-any_of("X")) %>%
  mutate(sample_vars = paste(Sample, CHROM, POS, REF, ALT, sep = "_"))

normalized_vcf_parsed_second = read.csv("/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/Tn5/second_round_20260130/results/variants_norm_WGS_20260219.csv") %>%
  select(-any_of("X")) %>%
  mutate(sample_vars = paste(Sample, CHROM, POS, REF, ALT, sep = "_"))

normalized_vcf_all = rbind(normalized_vcf_parsed, normalized_vcf_parsed_second)

designs_all = designs_all %>%
  mutate(
    edited = ifelse(sample_vars %in% normalized_vcf_all$sample_vars, 1, 0),
    assayed_colony = ifelse(!is.na(target_basecov) & target_basecov >= 2, 1, 0),
    edited_colony = ifelse(assayed_colony == 1 & edited == 1, 1, 0)
  )

colonies_with_coverage = designs_all %>% filter(target_basecov >= 2)

designs_assayed = designs %>%
  filter(uid %in% designs_with_target_cov) %>%
  mutate(
    edited = ifelse(uid %in% (colonies_with_coverage %>% filter(edited == 1))$uid, 1, 0),
    pam_dist = pmax(-3, pmin(pam_dist_var_start, pam_dist_var_end))
  )

designs_assayed = designs_assayed %>%
  mutate(pam_dist_bin = case_when(
    pam_dist <= 0 ~ "PAM",
    pam_dist <= 6 ~ "1-6",
    pam_dist <= 11 ~ "7-11",
    pam_dist <= 17 ~ "12-17",
    pam_dist <= 23 ~ "18-23",
    TRUE ~ ">23"
  ))

grouped_pam = designs_assayed %>%
  group_by(pam) %>%
  summarize(total_guides = n(), edited_guides = sum(edited), .groups = "drop") %>%
  mutate(edited_percentage = edited_guides / total_guides * 100)

grouped_pam_dist = designs_assayed %>%
  group_by(pam_dist) %>%
  summarize(total_guides = n(), edited_guides = sum(edited), .groups = "drop") %>%
  mutate(edited_percentage = edited_guides / total_guides * 100)

grouped_pam_dist_bin = designs_assayed %>%
  group_by(pam_dist_bin) %>%
  summarize(total_guides = n(), edited_guides = sum(edited), .groups = "drop") %>%
  mutate(
    edited_percentage = edited_guides / total_guides * 100,
    unedited = 100 - edited_percentage
  )

plot_pam_dist = grouped_pam_dist_bin %>%
  select(pam_dist_bin, edited_percentage, unedited) %>%
  pivot_longer(-pam_dist_bin, names_to = "status", values_to = "pct") %>%
  mutate(
    pam_dist_bin = factor(pam_dist_bin, levels = c("PAM", "1-6", "7-11", "12-17", "18-23")),
    status = recode(status, edited_percentage = "Edited", unedited = "Unedited"),
    status = factor(status, levels = c("Edited", "Unedited"))
  )

plot_pam_dist_svg = ggplot(plot_pam_dist, aes(pam_dist_bin, pct, fill = status)) +
  geom_col(width = 0.92, position = position_stack(reverse = TRUE)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_fill_manual(values = c(Edited = "#E69F00", Unedited = "#D9D9D9")) +
  labs(x = "PAM distance bin", y = "Editing (%)", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(hjust = 1),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1
  )

save_plot(plot_pam_dist_svg, file.path(plots_dir, "GW_pam_dist.svg"), width = 5, height = 5)

designs_all$delta_cov = designs_all$mean_cov - designs_all$target_basecov

coverage_qc_original = ggplot(designs_all %>% filter(mean_cov < 40), aes(x = mean_cov, y = target_basecov)) +
  geom_point(aes(color = as.factor(edited))) +
  xlim(0, 40) +
  ylim(0, 40) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Mean Coverage vs Target Base Coverage", x = "Target Base Coverage", y = "Mean Coverage") +
  theme_minimal() +
  geom_text(aes(label = ifelse(delta_cov > 10, colony_ID, "")), vjust = -1, position = position_jitter(width = 0.1, height = 0.1)) +
  theme(aspect.ratio = 1)

save_plot(coverage_qc_original, file.path(plots_dir, "coverage_qc_original_style.pdf"), width = 5.5, height = 5.5)

designs_assayed = designs_assayed %>%
  mutate(gc_content_bin = case_when(
    gc_content < 0.3 ~ "<30",
    gc_content < 0.4 ~ "30-40",
    gc_content < 0.5 ~ "40-50",
    gc_content < 0.6 ~ "50-60",
    gc_content < 0.7 ~ "60-70",
    TRUE ~ ">=70"
  ))

desings_assayed_gc_grouped = designs_assayed %>%
  group_by(gc_content_bin) %>%
  summarize(total_guides = n(), edited_guides = sum(edited), .groups = "drop") %>%
  mutate(edited_percentage = edited_guides / total_guides * 100)

designs_assayed_optimal = designs_assayed %>%
  filter(pam_dist_bin %in% c("1-6", "7-11", "12-17"), pam != "TTTG", max_A < 5)

sum(designs_assayed_optimal$edited) / nrow(designs_assayed_optimal)

plot_maxA_original = ggplot(designs_assayed, aes(x = max_A, y = edited)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Edited Status by Max A", x = "Max A", y = "Edited Status") +
  theme_minimal()

save_plot(plot_maxA_original, file.path(plots_dir, "maxA_original_style.pdf"), width = 5, height = 5)

grouped_maxA = designs_assayed %>%
  group_by(max_A) %>%
  summarize(total_guides = n(), edited_guides = sum(edited), .groups = "drop") %>%
  mutate(edited_percentage = edited_guides / total_guides * 100)

grouped_maxT = designs_assayed %>%
  group_by(max_T) %>%
  summarize(total_guides = n(), edited_guides = sum(edited), .groups = "drop") %>%
  mutate(edited_percentage = edited_guides / total_guides * 100)

model_optimal = glm(edited ~ pam + pam_dist_bin + max_A + max_T, family = binomial, data = designs_assayed)
summary(model_optimal)

designs_assayed_optimal = designs_assayed %>%
  filter(!pam_dist_bin %in% c("PAM", "18-23"), pam != "TTTG", max_A < 4)

## DeepCpf1
ref = readDNAStringSet("/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/Tn5/Analysis/yKR61_fastq/saccharomyces_cerevisiae_sequence.fasta")
names(ref) <- sub("\\s.*", "", names(ref))

extract_guide_context = function(chr, pos, strand, guide_seq, ref, flank_up = 3L, flank_down = 1L, search_pad = 80L) {
  chr_seq = ref[[as.character(chr)]]
  chr_len = length(chr_seq)

  win_start = max(1L, pos - search_pad)
  win_end = min(chr_len, pos + search_pad + nchar(guide_seq))
  local_seq = subseq(chr_seq, start = win_start, end = win_end)

  query = if (strand == "+") {
    DNAString(guide_seq)
  } else {
    reverseComplement(DNAString(guide_seq))
  }

  hits = matchPattern(query, local_seq, fixed = TRUE)
  if (length(hits) == 0) return(NA_character_)

  hit_start = win_start + start(hits)[1] - 1L
  hit_end = win_start + end(hits)[1] - 1L

  if (strand == "+") {
    out = subseq(chr_seq, start = max(1L, hit_start - flank_up), end = min(chr_len, hit_end + flank_down))
    as.character(out)
  } else {
    out = subseq(chr_seq, start = max(1L, hit_start - flank_down), end = min(chr_len, hit_end + flank_up))
    as.character(reverseComplement(out))
  }
}

designs_out_deep_model = designs %>%
  mutate(
    pos = as.integer(str_match(guide_id, "^pos_0*([0-9]+)_")[, 2]),
    strand = str_match(guide_id, "_([+-])_strand$")[, 2],
    target_seq = pmap_chr(list(CHR, pos, strand, guide_seq), ~ extract_guide_context(..1, ..2, ..3, ..4, ref)),
    guide_core = substr(target_seq, 4, 30),
    matches_guide = guide_core == guide_seq
  )

deep_cpf_input = designs_out_deep_model %>% select(guide_id, target_seq)
write.csv(deep_cpf_input, file.path(inputs_dir, "input_deepcpf1_20260319.csv"), row.names = FALSE)

deep_cpf_results = read.csv(file.path(project_root, "input_deepcpf1_20260319.scored.csv"))
designs_assayed_deepcpf = designs_assayed %>% left_join(deep_cpf_results, by = "guide_id")

check_pred_vs_observed = function(data, pred_col = "y_pred", obs_col = "edited") {
  df = data %>%
    dplyr::select(pred = dplyr::all_of(pred_col), obs = dplyr::all_of(obs_col)) %>%
    dplyr::filter(!is.na(pred), !is.na(obs)) %>%
    dplyr::mutate(obs = dplyr::case_when(
      obs %in% c(TRUE, "TRUE", "Yes", "YES", "yes", 1, "1") ~ 1,
      obs %in% c(FALSE, "FALSE", "No", "NO", "no", 0, "0") ~ 0,
      TRUE ~ NA_real_
    )) %>%
    dplyr::filter(!is.na(obs))

  if (nrow(df) == 0) stop("No valid rows after filtering.")
  if (length(unique(df$obs)) < 2) stop("Observed column must contain both edited and non-edited values.")

  pearson_test = suppressWarnings(cor.test(df$pred, df$obs, method = "pearson"))
  spearman_test = suppressWarnings(cor.test(df$pred, df$obs, method = "spearman"))
  auc_val = as.numeric(pROC::roc(response = df$obs, predictor = df$pred, quiet = TRUE)$auc)
  glm_fit = glm(obs ~ pred, data = df, family = binomial())

  tibble(
    n = nrow(df),
    edited_n = sum(df$obs == 1),
    unedited_n = sum(df$obs == 0),
    pearson_r = unname(pearson_test$estimate),
    pearson_p = pearson_test$p.value,
    spearman_rho = unname(spearman_test$estimate),
    spearman_p = spearman_test$p.value,
    auc = auc_val,
    glm_beta = unname(coef(glm_fit)[["pred"]]),
    glm_odds_ratio = exp(unname(coef(glm_fit)[["pred"]])),
    glm_p = summary(glm_fit)$coefficients["pred", "Pr(>|z|)"]
  )
}

deepcpf1_prediction_summary = check_pred_vs_observed(designs_assayed_deepcpf, pred_col = "y_pred", obs_col = "edited")

## DR folding
DR_seq = "TAATTTCTACTCTTGTAGAT"
designs_vienna = designs %>%
  select(uid, guide_nopam) %>%
  mutate(DR_spacer = paste0(DR_seq, guide_nopam)) %>%
  select(uid, DR_spacer)

write.csv(designs_vienna, file.path(inputs_dir, "guide_table_rna_20260322.csv"), row.names = FALSE)

rna_fold_guides = read.csv(
  "/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/Tn5/Analysis/annotations_guides/guides_dr_scores_20260325.tsv",
  header = TRUE,
  sep = "\t"
)

designs_assayed = left_join(designs_assayed, rna_fold_guides, by = c("uid" = "guide_id"))
designs_assayed_deepcpf = left_join(designs_assayed_deepcpf, rna_fold_guides, by = c("uid" = "guide_id"))

designs_assayed_unpaired = designs_assayed %>%
  filter(structure == ".....(((((....)))))........................")

glm(edited ~ mfe_kcal_mol * pam, data = designs_assayed, family = binomial) %>% summary

designs_assayed = designs_assayed %>%
  mutate(
    gc_content_20nt = (str_count(substr(guide_nopam, 1, 20), "G") + str_count(substr(guide_nopam, 1, 20), "C")) / 20,
    gc_content_20nt_bin = case_when(
      gc_content_20nt < 0.3 ~ "<30%",
      gc_content_20nt < 0.4 ~ "30-39%",
      gc_content_20nt < 0.5 ~ "40-49%",
      gc_content_20nt < 0.6 ~ "50-59%",
      TRUE ~ ">=60%"
    ),
    gc_content_18nt = (str_count(substr(guide_nopam, 1, 18), "G") + str_count(substr(guide_nopam, 1, 18), "C")) / 18,
    gc_content_18nt_bin = case_when(
      gc_content_18nt < 0.3 ~ "<30%",
      gc_content_18nt < 0.4 ~ "30-39%",
      gc_content_18nt < 0.5 ~ "40-49%",
      gc_content_18nt < 0.6 ~ "50-59%",
      TRUE ~ ">=60%"
    ),
    donor_wt_gc = gc_fraction(donor_wt),
    donor_wt_gc_bin = gc_fraction_bin(donor_wt_gc),
    gc_content_bin = gc_fraction_bin(gc_content),
    dr_perturbed = tolower(as.character(dr_perturbed)) %in% c("yes", "true", "1"),
    dr_canonical_structure = structure == ".....(((((....)))))........................",
    dr_spacer_pair_any = dr_nt_paired_to_spacer > 0
  )

grouped_gc_content_20nt_bin = designs_assayed %>%
  group_by(gc_content_20nt_bin) %>%
  summarize(total_guides = n(), edited_guides = sum(edited), .groups = "drop") %>%
  mutate(edited_percentage = edited_guides / total_guides * 100)

grouped_gc_content_18nt_bin = designs_assayed %>%
  group_by(gc_content_18nt_bin) %>%
  summarize(total_guides = n(), edited_guides = sum(edited), .groups = "drop") %>%
  mutate(edited_percentage = edited_guides / total_guides * 100)

designs_assayed_deepcpf = designs_assayed_deepcpf %>%
  mutate(
    ntt_pattern = substr(target_seq, 3, 7),
    ntt_pattern_upstream = substr(target_seq, 3, 3)
  )

designs_assayed_deepcpf$ntt_pattern_upstream = factor(designs_assayed_deepcpf$ntt_pattern_upstream, levels = c("C", "A", "G", "T"))

grouped_ntt_pattern = designs_assayed_deepcpf %>%
  group_by(ntt_pattern) %>%
  summarize(total_guides = n(), edited_guides = sum(edited), .groups = "drop") %>%
  mutate(edited_percentage = edited_guides / total_guides * 100)

grouped_ntt_pattern_upstream = designs_assayed_deepcpf %>%
  group_by(ntt_pattern_upstream) %>%
  summarize(total_guides = n(), edited_guides = sum(edited), .groups = "drop") %>%
  mutate(edited_percentage = edited_guides / total_guides * 100)

fit_upstream = glm(edited ~ ntt_pattern_upstream, data = designs_assayed_deepcpf, family = binomial)

designs_assayed_optimal_ntt = designs_assayed %>%
  left_join(designs_assayed_deepcpf %>% select(uid, ntt_pattern, ntt_pattern_upstream), by = "uid") %>%
  filter(!pam_dist_bin %in% c("PAM", "18-23"), pam != "TTTG", max_A < 4, ntt_pattern_upstream == "T")

## design-level concordance in the original style
colonies_multiple_barcodes = colonies_with_coverage %>%
  group_by(uid) %>%
  mutate(num_barcodes = n(), fraction_edited = sum(edited) / n())

concordance_summary = colonies_multiple_barcodes %>%
  group_by(uid) %>%
  summarize(
    num_barcodes = first(num_barcodes),
    fraction_edited = first(fraction_edited),
    concordance_category = case_when(
      fraction_edited == 0 ~ "All unedited",
      fraction_edited == 1 ~ "All edited",
      TRUE ~ "Mixed"
    ),
    .groups = "drop"
  )

## notebook-style design-rule tables
score_sv_sites = read.delim("/Users/u0174312/Downloads/s288c_NGG_PAM_sites_editing_SCORE.tsv", sep = "\t", header = TRUE) %>%
  transmute(
    chr = paste0("chr", chr),
    pos = as.integer(pos),
    score_sv = as.numeric(`SCORE.SV.minmax`),
    guide = guide,
    id = id
  )

extra_design_flags = design_annotations_extra %>%
  select(uid, in_sv_hotspot, exact_high_score, hotspot_high_score) %>%
  distinct()

designs_all_annotated = designs_all %>%
  left_join(select(deep_cpf_results, guide_id, target_seq, y_pred), by = "guide_id") %>%
  left_join(rna_fold_guides, by = c("uid" = "guide_id")) %>%
  left_join(extra_design_flags, by = "uid") %>%
  left_join(
    variant_annotations %>% select(chr, pos, ref, alt, region_class, snpeff_primary_impact, snpeff_primary_effect, spcas9_targetable),
    by = c("CHR" = "chr", "POS" = "pos", "REF" = "ref", "ALT" = "alt")
  ) %>%
  annotate_nearest_score_sv("CHR", "POS", score_sv_sites) %>%
  mutate(
    pam_dist = pmax(-3, pmin(pam_dist_var_start, pam_dist_var_end)),
    pam_dist_bin = case_when(
      pam_dist <= 0 ~ "PAM",
      pam_dist <= 6 ~ "1-6",
      pam_dist <= 11 ~ "7-11",
      pam_dist <= 17 ~ "12-17",
      pam_dist <= 23 ~ "18-23",
      TRUE ~ ">23"
    ),
    pam_dist_bin = factor(pam_dist_bin, levels = c("1-6", "7-11", "12-17", "18-23", "PAM", ">23")),
    variant_class = ifelse(vardiff == 0, "SNV", "Indel"),
    ntt_pattern = substr(target_seq, 3, 7),
    ntt_pattern_upstream = substr(target_seq, 3, 3),
    full_pam_5mer = paste0(ntt_pattern_upstream, pam),
    guide_gc_23 = gc_content,
    guide_gc_20 = gc_fraction(substr(guide_nopam, 1, 20)),
    guide_gc_18 = gc_fraction(substr(guide_nopam, 1, 18)),
    donor_wt_gc = gc_fraction(donor_wt),
    guide_gc_23_bin = gc_fraction_bin(guide_gc_23),
    guide_gc_20_bin = gc_fraction_bin(guide_gc_20),
    guide_gc_18_bin = gc_fraction_bin(guide_gc_18),
    donor_wt_gc_bin = gc_fraction_bin(donor_wt_gc),
    dr_perturbed = tolower(as.character(dr_perturbed)) %in% c("yes", "true", "1"),
    dr_canonical_structure = structure == ".....(((((....)))))........................",
    dr_spacer_pair_any = dr_nt_paired_to_spacer > 0
  )

deepcpf1_quartiles = designs_all_annotated %>%
  distinct(uid, y_pred) %>%
  mutate(
    deepcpf1_quartile = ntile(y_pred, 4),
    deepcpf1_quartile = factor(deepcpf1_quartile, levels = 1:4, labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest")),
    high_deepcpf1 = y_pred >= median(y_pred, na.rm = TRUE)
  )

designs_all_annotated = designs_all_annotated %>%
  left_join(deepcpf1_quartiles, by = c("uid", "y_pred")) %>%
  mutate(
    preferred_pam = pam %in% c("TTTA", "TTTC"),
    preferred_distance = pam_dist_bin %in% c("1-6", "7-11", "12-17"),
    preferred_rule_count = as.integer(preferred_pam) + as.integer(preferred_distance) + as.integer(high_deepcpf1),
    recommended_design = preferred_rule_count == 3
  )

design_outcomes = designs_all_annotated %>%
  group_by(uid) %>%
  summarize(
    round = first(round),
    guide_id = first(guide_id),
    pam = first(pam),
    pam_dist = first(pam_dist),
    pam_dist_bin = first(pam_dist_bin),
    deepcpf1_score = first(y_pred),
    deepcpf1_quartile = first(deepcpf1_quartile),
    pam_upstream_nt = first(ntt_pattern_upstream),
    full_pam_5mer = first(full_pam_5mer),
    ntttv_pattern = first(ntt_pattern),
    gc_content = first(gc_content),
    guide_gc_20 = first(guide_gc_20),
    guide_gc_18 = first(guide_gc_18),
    donor_wt_gc = first(donor_wt_gc),
    guide_gc_23_bin = first(guide_gc_23_bin),
    guide_gc_20_bin = first(guide_gc_20_bin),
    guide_gc_18_bin = first(guide_gc_18_bin),
    donor_wt_gc_bin = first(donor_wt_gc_bin),
    max_A = first(max_A),
    max_T = first(max_T),
    num_T = first(num_T),
    vardiff = first(vardiff),
    variant_class = first(variant_class),
    region_class = first(region_class),
    snpeff_primary_impact = first(snpeff_primary_impact),
    snpeff_primary_effect = first(snpeff_primary_effect),
    spcas9_targetable = first(spcas9_targetable),
    in_sv_hotspot = first(in_sv_hotspot),
    exact_high_score = first(exact_high_score),
    hotspot_high_score = first(hotspot_high_score),
    nearest_score_sv = first(nearest_score_sv),
    nearest_score_sv_pos = first(nearest_score_sv_pos),
    nearest_score_sv_guide = first(nearest_score_sv_guide),
    dr_mfe_kcal_mol = first(mfe_kcal_mol),
    dr_nt_paired_to_dr = first(dr_nt_paired_to_dr),
    dr_nt_paired_to_spacer = first(dr_nt_paired_to_spacer),
    dr_nt_unpaired = first(dr_nt_unpaired),
    dr_perturbed = first(dr_perturbed),
    dr_score = first(score),
    dr_canonical_structure = first(dr_canonical_structure),
    dr_spacer_pair_any = first(dr_spacer_pair_any),
    preferred_pam = first(preferred_pam),
    preferred_distance = first(preferred_distance),
    high_deepcpf1 = first(high_deepcpf1),
    preferred_rule_count = first(preferred_rule_count),
    recommended_design = first(recommended_design),
    picked_colonies = n(),
    assayed_colonies = sum(assayed_colony, na.rm = TRUE),
    edited_colonies = sum(edited_colony, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    assayed_design = assayed_colonies > 0,
    edited_design = edited_colonies > 0,
    colony_edit_rate = ifelse(assayed_colonies > 0, edited_colonies / assayed_colonies, NA_real_)
  )

design_coordinates = designs_all_annotated %>%
  group_by(uid) %>%
  summarize(
    chr = first(as.character(CHR)),
    pos = first(as.integer(POS)),
    .groups = "drop"
  )

chromosome_order = names(ref)[names(ref) != "chrmt" & names(ref) %in% unique(design_coordinates$chr)]
chromosome_gap_bp = 50000L

chromosome_layout = tibble(
  chr = chromosome_order,
  chr_len = as.integer(width(ref[chromosome_order]))
) %>%
  mutate(
    chr_gap = chromosome_gap_bp,
    genome_start = lag(cumsum(chr_len + chr_gap), default = 0),
    genome_end = genome_start + chr_len
  )

chromosome_total_span = max(chromosome_layout$genome_end) + chromosome_gap_bp

chromosome_layout = chromosome_layout %>%
  mutate(
    genome_mid = (genome_start + genome_end) / 2,
    chr_label = sub("^chr", "", chr),
    chromosome_label_angle = 90 - 360 * genome_mid / chromosome_total_span,
    chromosome_label_hjust = ifelse(chromosome_label_angle < -90, 1, 0),
    chromosome_label_angle = ifelse(chromosome_label_angle < -90, chromosome_label_angle + 180, chromosome_label_angle)
  )

chromosome_bin_width = 100000L

chromosome_bin_layout = map_dfr(
  chromosome_order,
  function(chr_name) {
    chr_len = chromosome_layout$chr_len[chromosome_layout$chr == chr_name]
    bin_starts = seq(1L, chr_len, by = chromosome_bin_width)

    tibble(
      chr = chr_name,
      bin_id = seq_along(bin_starts),
      bin_start = bin_starts,
      bin_end = pmin(bin_starts + chromosome_bin_width - 1L, chr_len)
    )
  }
) %>%
  left_join(chromosome_layout, by = "chr") %>%
  mutate(
    bin_width_bp = bin_end - bin_start + 1L,
    genome_bin_start = genome_start + bin_start - 1L,
    genome_bin_end = genome_start + bin_end,
    bin_plot_width_bp = pmin(
      pmax(as.integer(round(bin_width_bp * 0.72)), 1500L),
      pmax(bin_width_bp - 2L, 1L)
    ),
    bin_plot_start = genome_bin_start + floor((bin_width_bp - bin_plot_width_bp) / 2),
    bin_plot_end = pmin(genome_bin_end - 1L, bin_plot_start + bin_plot_width_bp)
  )

chromosome_circos_summary = design_outcomes %>%
  left_join(design_coordinates, by = "uid") %>%
  filter(assayed_design, chr %in% chromosome_order) %>%
  mutate(
    bin_id = floor((pos - 1L) / chromosome_bin_width) + 1L
  ) %>%
  group_by(chr, bin_id) %>%
  summarize(
    assayed_designs = n(),
    edited_designs = sum(edited_design),
    .groups = "drop"
  ) %>%
  right_join(chromosome_bin_layout, by = c("chr", "bin_id")) %>%
  arrange(factor(chr, levels = chromosome_order), bin_id) %>%
  mutate(
    assayed_designs = replace_na(assayed_designs, 0L),
    edited_designs = replace_na(edited_designs, 0L)
  ) %>%
  bind_cols(
    binom_wilson(.$edited_designs, pmax(.$assayed_designs, 1L)) %>%
      dplyr::rename(
        design_success_rate = rate,
        design_success_lower = lower,
        design_success_upper = upper
      )
  ) %>%
  mutate(
    design_success_rate = ifelse(assayed_designs > 0, design_success_rate, NA_real_),
    design_success_lower = ifelse(assayed_designs > 0, design_success_lower, NA_real_),
    design_success_upper = ifelse(assayed_designs > 0, design_success_upper, NA_real_),
    assayed_count_scaled = assayed_designs / max(assayed_designs)
  )

stage_summary = bind_rows(
  tibble(
    round = "overall",
    picked_designs = nrow(design_outcomes),
    assayed_designs = sum(design_outcomes$assayed_design),
    edited_designs = sum(design_outcomes$edited_design),
    picked_colonies = sum(design_outcomes$picked_colonies),
    assayed_colonies = sum(design_outcomes$assayed_colonies),
    edited_colonies = sum(design_outcomes$edited_colonies)
  ),
  design_outcomes %>%
    group_by(round) %>%
    summarize(
      picked_designs = n(),
      assayed_designs = sum(assayed_design),
      edited_designs = sum(edited_design),
      picked_colonies = sum(picked_colonies),
      assayed_colonies = sum(assayed_colonies),
      edited_colonies = sum(edited_colonies),
      .groups = "drop"
    )
)

pam_summary = summarize_edit_rates(design_outcomes, "pam") %>%
  mutate(pam = factor(pam, levels = c("TTTC", "TTTA", "TTTG")))

pam_distance_summary = summarize_edit_rates(design_outcomes, "pam_dist_bin") %>%
  mutate(pam_dist_bin = factor(pam_dist_bin, levels = c("1-6", "7-11", "12-17", "18-23", "PAM")))

pam_heatmap_summary = summarize_edit_rates(design_outcomes, c("pam", "pam_dist_bin")) %>%
  filter(!is.na(pam_dist_bin), pam_dist_bin != ">23")

deepcpf1_summary = summarize_edit_rates(design_outcomes, "deepcpf1_quartile")
pam_upstream_summary = summarize_edit_rates(design_outcomes, "pam_upstream_nt") %>%
  mutate(pam_upstream_nt = factor(pam_upstream_nt, levels = c("C", "A", "G", "T")))
full_pam_5mer_summary = summarize_edit_rates(design_outcomes, "full_pam_5mer") %>% arrange(full_pam_5mer)
ntttv_pattern_summary = summarize_edit_rates(design_outcomes, "ntttv_pattern") %>% arrange(desc(total_guides), ntttv_pattern)
rule_count_summary = summarize_edit_rates(design_outcomes, "preferred_rule_count") %>%
  mutate(preferred_rule_count = factor(preferred_rule_count, levels = 0:3))

dr_perturbed_summary = design_outcomes %>%
  mutate(
    dr_perturbed_flag = ifelse(dr_perturbed, "Perturbed", "Not perturbed"),
    dr_perturbed_flag = factor(dr_perturbed_flag, levels = c("Not perturbed", "Perturbed"))
  ) %>%
  summarize_edit_rates("dr_perturbed_flag")

dr_canonical_summary = design_outcomes %>%
  mutate(
    dr_canonical_flag = ifelse(dr_canonical_structure, "Canonical", "Non-canonical"),
    dr_canonical_flag = factor(dr_canonical_flag, levels = c("Canonical", "Non-canonical"))
  ) %>%
  summarize_edit_rates("dr_canonical_flag")

dr_spacer_pair_summary = design_outcomes %>%
  mutate(
    dr_spacer_pair_flag = ifelse(dr_spacer_pair_any, "Any DR-spacer pairing", "No DR-spacer pairing"),
    dr_spacer_pair_flag = factor(dr_spacer_pair_flag, levels = c("No DR-spacer pairing", "Any DR-spacer pairing"))
  ) %>%
  summarize_edit_rates("dr_spacer_pair_flag")

guide_gc_23_summary = summarize_edit_rates(design_outcomes, "guide_gc_23_bin") %>%
  mutate(guide_gc_23_bin = factor(guide_gc_23_bin, levels = c("<30%", "30-39%", "40-49%", "50-59%", ">=60%")))

guide_gc_20_summary = summarize_edit_rates(design_outcomes, "guide_gc_20_bin") %>%
  mutate(guide_gc_20_bin = factor(guide_gc_20_bin, levels = c("<30%", "30-39%", "40-49%", "50-59%", ">=60%")))

guide_gc_18_summary = summarize_edit_rates(design_outcomes, "guide_gc_18_bin") %>%
  mutate(guide_gc_18_bin = factor(guide_gc_18_bin, levels = c("<30%", "30-39%", "40-49%", "50-59%", ">=60%")))

donor_wt_gc_summary = summarize_edit_rates(design_outcomes, "donor_wt_gc_bin") %>%
  mutate(donor_wt_gc_bin = factor(donor_wt_gc_bin, levels = c("<30%", "30-39%", "40-49%", "50-59%", ">=60%")))

max_A_summary = summarize_edit_rates(design_outcomes, "max_A") %>%
  mutate(max_A = factor(max_A, levels = sort(unique(max_A))))

max_T_summary = summarize_edit_rates(design_outcomes, "max_T") %>%
  mutate(max_T = factor(max_T, levels = sort(unique(max_T))))

region_class_summary = design_outcomes %>%
  mutate(region_class = replace_na(region_class, "unknown")) %>%
  summarize_edit_rates("region_class") %>%
  mutate(region_class = factor(region_class, levels = c("coding", "promoter", "intergenic", "unknown")))

hotspot_summary = design_outcomes %>%
  mutate(
    hotspot_flag = ifelse(coalesce(in_sv_hotspot, FALSE), "SV hotspot", "Not hotspot"),
    hotspot_flag = factor(hotspot_flag, levels = c("Not hotspot", "SV hotspot"))
  ) %>%
  summarize_edit_rates("hotspot_flag")

multi_colony_consistency_summary = design_outcomes %>%
  filter(assayed_design, assayed_colonies >= 2) %>%
  mutate(
    consistency = case_when(
      edited_colonies == 0 ~ "None edited",
      edited_colonies == assayed_colonies ~ "All edited",
      TRUE ~ "Mixed"
    ),
    consistency = factor(consistency, levels = c("All edited", "Mixed", "None edited"))
  ) %>%
  count(assayed_colonies, consistency, name = "design_count") %>%
  group_by(assayed_colonies) %>%
  mutate(total_designs = sum(design_count), fraction = design_count / total_designs) %>%
  ungroup()

design_rule_summary = bind_rows(
  summarize_rule_comparison(
    design_outcomes %>% filter(assayed_design, preferred_pam),
    design_outcomes %>% filter(assayed_design, !preferred_pam),
    "Primary",
    "Prefer TTTA or TTTC PAMs over TTTG",
    "TTTA or TTTC",
    "TTTG",
    "Consistent in both marginal summaries and the unified binary cassette-level model."
  ),
  summarize_rule_comparison(
    design_outcomes %>% filter(assayed_design, preferred_distance),
    design_outcomes %>% filter(assayed_design, !preferred_distance),
    "Primary",
    "Place the edit 1-17 nt from the PAM",
    "1-17 nt from PAM",
    "PAM-overlapping or 18-23 nt",
    "PAM-overlapping and 18-23 nt were the main distance failures."
  ),
  summarize_rule_comparison(
    design_outcomes %>% filter(assayed_design, deepcpf1_quartile == "Q4 highest"),
    design_outcomes %>% filter(assayed_design, deepcpf1_quartile == "Q1 lowest"),
    "Primary",
    "Prioritize guides with higher DeepCpf1 scores",
    "DeepCpf1 Q4",
    "DeepCpf1 Q1",
    "Treat the score as a prioritization signal, not an absolute cutoff."
  ),
  summarize_rule_comparison(
    design_outcomes %>% filter(assayed_design, recommended_design),
    design_outcomes %>% filter(assayed_design, !recommended_design),
    "Primary",
    "Best simple bundle: preferred PAM + 1-17 nt distance + DeepCpf1 above the median",
    "All 3 preferred rules",
    "Anything else",
    "This is the cleanest simple design bundle in the dataset."
  ),
  summarize_rule_comparison(
    design_outcomes %>% filter(assayed_design, max_A < 5),
    design_outcomes %>% filter(assayed_design, max_A >= 5),
    "Secondary",
    "Use long A-runs only as a tie-breaker filter",
    "max_A < 5",
    "max_A >= 5",
    "The drop for long A-runs is visible, but it is not a primary rule."
  ),
  summarize_rule_comparison(
    design_outcomes %>% filter(assayed_design, !coalesce(in_sv_hotspot, FALSE)),
    design_outcomes %>% filter(assayed_design, coalesce(in_sv_hotspot, FALSE)),
    "Context",
    "Avoid SV hotspot loci if you have equivalent alternatives",
    "Not in SV hotspot",
    "SV hotspot",
    "Useful as a contextual filter, not a primary cassette rule."
  )
)

variant_calls_all = normalized_vcf_all %>%
  transmute(
    sample = Sample,
    CHR = CHROM,
    POS = as.integer(POS),
    REF = REF,
    ALT = ALT,
    QUAL = as.numeric(QUAL),
    DP = as.numeric(DP),
    AF = as.numeric(frc_alt),
    TYPE = TYPE,
    sample_vars = paste(Sample, CHROM, POS, REF, ALT, sep = "_")
  )

intended_sites = designs_all_annotated %>%
  transmute(
    sample = sample,
    intended_sample_vars = sample_vars,
    assayed_colony = assayed_colony == 1,
    edited_colony = edited_colony == 1
  ) %>%
  distinct()

background_variants = variant_calls_all %>%
  left_join(intended_sites %>% select(sample, intended_sample_vars), by = "sample") %>%
  mutate(is_intended_variant = sample_vars == intended_sample_vars) %>%
  select(-intended_sample_vars) %>%
  annotate_nearest_score_sv("CHR", "POS", score_sv_sites) %>%
  mutate(
    high_conf_non_target = !is_intended_variant & QUAL >= 20 & AF >= 50 & DP >= 4,
    high_score_sv_site = nearest_score_sv >= 0.9,
    variant_id = paste(CHR, POS, REF, ALT, sep = ":")
  )

high_conf_non_target_variants = background_variants %>% filter(high_conf_non_target)

background_colony_summary = intended_sites %>%
  left_join(
    high_conf_non_target_variants %>%
      group_by(sample) %>%
      summarize(
        high_conf_non_target_count = n(),
        max_non_target_score_sv = max(nearest_score_sv, na.rm = TRUE),
        mean_non_target_score_sv = mean(nearest_score_sv, na.rm = TRUE),
        frac_non_target_high_score = mean(high_score_sv_site, na.rm = TRUE),
        max_non_target_qual = max(QUAL, na.rm = TRUE),
        max_non_target_af = max(AF, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "sample"
  ) %>%
  mutate(high_conf_non_target_count = coalesce(high_conf_non_target_count, 0L))

background_edit_status_summary = background_colony_summary %>%
  mutate(target_edit_status = ifelse(edited_colony, "Edited target", "Unedited target")) %>%
  group_by(target_edit_status) %>%
  summarize(
    colony_count = n(),
    mean_high_conf_non_target_count = mean(high_conf_non_target_count, na.rm = TRUE),
    median_high_conf_non_target_count = median(high_conf_non_target_count, na.rm = TRUE),
    frac_with_any_high_conf_non_target = mean(high_conf_non_target_count > 0, na.rm = TRUE),
    mean_max_non_target_score_sv = mean(max_non_target_score_sv, na.rm = TRUE),
    .groups = "drop"
  )

recurrent_background_variants = high_conf_non_target_variants %>%
  count(variant_id, CHR, POS, REF, ALT, nearest_score_sv, name = "sample_count", sort = TRUE) %>%
  filter(sample_count >= 2)

design_background_summary = designs_all_annotated %>%
  filter(assayed_colony == 1) %>%
  distinct(uid, sample) %>%
  left_join(
    background_colony_summary %>% select(sample, high_conf_non_target_count),
    by = "sample"
  ) %>%
  mutate(high_conf_non_target_count = coalesce(high_conf_non_target_count, 0L)) %>%
  group_by(uid) %>%
  summarize(
    other_variant_design = any(high_conf_non_target_count > 0),
    assayed_colonies_with_other_variants = sum(high_conf_non_target_count > 0),
    .groups = "drop"
  )

design_outcomes = design_outcomes %>%
  left_join(design_background_summary, by = "uid") %>%
  mutate(
    other_variant_design = coalesce(other_variant_design, FALSE),
    assayed_colonies_with_other_variants = coalesce(assayed_colonies_with_other_variants, 0L),
    outcome_state = case_when(
      other_variant_design ~ "Other variants observed",
      edited_design ~ "Intended edit only",
      TRUE ~ "Unedited"
    ),
    outcome_state = factor(
      outcome_state,
      levels = c("Intended edit only", "Other variants observed", "Unedited")
    )
  )

stacked_parameter_data = design_outcomes %>%
  mutate(
    pam = factor(pam, levels = c("TTTC", "TTTA", "TTTG")),
    pam_upstream_nt = factor(pam_upstream_nt, levels = c("C", "A", "G", "T")),
    pam_dist_bin = factor(pam_dist_bin, levels = c("PAM", "1-6", "7-11", "12-17", "18-23", ">23")),
    deepcpf1_quartile = factor(deepcpf1_quartile, levels = c("Q1 lowest", "Q2", "Q3", "Q4 highest")),
    guide_gc_23_bin = factor(guide_gc_23_bin, levels = c("<30%", "30-39%", "40-49%", "50-59%", ">=60%")),
    guide_gc_20_bin = factor(guide_gc_20_bin, levels = c("<30%", "30-39%", "40-49%", "50-59%", ">=60%")),
    guide_gc_18_bin = factor(guide_gc_18_bin, levels = c("<30%", "30-39%", "40-49%", "50-59%", ">=60%")),
    donor_wt_gc_bin = factor(donor_wt_gc_bin, levels = c("<30%", "30-39%", "40-49%", "50-59%", ">=60%")),
    dr_perturbed_flag = factor(ifelse(dr_perturbed, "Perturbed", "Not perturbed"), levels = c("Not perturbed", "Perturbed")),
    dr_canonical_flag = factor(ifelse(dr_canonical_structure, "Canonical", "Non-canonical"), levels = c("Canonical", "Non-canonical")),
    dr_spacer_pair_flag = factor(ifelse(dr_spacer_pair_any, "Any DR-spacer pairing", "No DR-spacer pairing"), levels = c("No DR-spacer pairing", "Any DR-spacer pairing"))
  )

stacked_parameter_specs = tribble(
  ~parameter_col, ~parameter_pretty, ~parameter_group,
  "deepcpf1_quartile", "DeepCpf1 quartile", "Core cassette parameters",
  "pam", "PAM", "Core cassette parameters",
  "pam_upstream_nt", "5' PAM base", "Core cassette parameters",
  "pam_dist_bin", "PAM distance bin", "Core cassette parameters",
  "guide_gc_23_bin", "Guide GC (23 nt)", "GC and DR parameters",
  "guide_gc_20_bin", "Guide GC (first 20 nt)", "GC and DR parameters",
  "guide_gc_18_bin", "Guide GC (first 18 nt)", "GC and DR parameters",
  "donor_wt_gc_bin", "WT donor GC", "GC and DR parameters",
  "dr_perturbed_flag", "DR perturbation", "GC and DR parameters",
  "dr_canonical_flag", "DR structure class", "GC and DR parameters",
  "dr_spacer_pair_flag", "DR-spacer pairing", "GC and DR parameters"
)

stacked_design_outcome_summary = pmap_dfr(
  stacked_parameter_specs,
  function(parameter_col, parameter_pretty, parameter_group) {
    summarize_stacked_design_outcomes(stacked_parameter_data, parameter_col, parameter_pretty, parameter_group)
  }
)

unified_model_data = design_outcomes %>%
  filter(assayed_design) %>%
  mutate(
    round = factor(round, levels = c("round_1", "round_2")),
    pam = factor(pam, levels = c("TTTA", "TTTC", "TTTG")),
    pam_dist_bin = factor(pam_dist_bin, levels = c("1-6", "7-11", "12-17", "18-23", "PAM", ">23")),
    pam_upstream_nt = factor(pam_upstream_nt, levels = c("C", "A", "G", "T")),
    dr_spacer_pair_any = factor(dr_spacer_pair_any, levels = c(FALSE, TRUE))
  )

unified_model_main = glm(
  edited_design ~ pam + pam_upstream_nt + pam_dist_bin + deepcpf1_score + guide_gc_20 + donor_wt_gc + dr_spacer_pair_any,
  family = binomial(),
  data = unified_model_data
)

unified_model_pam_by_distance = glm(
  edited_design ~ pam + pam_upstream_nt + pam_dist_bin + deepcpf1_score + guide_gc_20 + donor_wt_gc + dr_spacer_pair_any + pam:pam_dist_bin,
  family = binomial(),
  data = unified_model_data
)

unified_model_score_by_distance = glm(
  edited_design ~ pam + pam_upstream_nt + pam_dist_bin + deepcpf1_score + guide_gc_20 + donor_wt_gc + dr_spacer_pair_any + deepcpf1_score:pam_dist_bin,
  family = binomial(),
  data = unified_model_data
)

unified_model_pam5_by_pam = glm(
  edited_design ~ pam + pam_upstream_nt + pam_dist_bin + deepcpf1_score + guide_gc_20 + donor_wt_gc + dr_spacer_pair_any + pam_upstream_nt:pam,
  family = binomial(),
  data = unified_model_data
)

unified_model_candidates = list(
  main_effects = unified_model_main,
  pam_by_distance = unified_model_pam_by_distance,
  score_by_distance = unified_model_score_by_distance,
  pam5_by_pam = unified_model_pam5_by_pam
)

interaction_model_summary = tibble(
  model = names(unified_model_candidates),
  formula = map_chr(unified_model_candidates, ~ paste(deparse(formula(.x)), collapse = "")),
  aic = map_dbl(unified_model_candidates, AIC),
  bic = map_dbl(unified_model_candidates, BIC),
  logLik = map_dbl(unified_model_candidates, ~ as.numeric(logLik(.x))),
  df = map_dbl(unified_model_candidates, ~ attr(logLik(.x), "df"))
) %>%
  mutate(
    delta_aic = aic - min(aic),
    delta_bic = bic - min(bic),
    lrt_p_vs_main = c(
      NA_real_,
      anova(unified_model_main, unified_model_pam_by_distance, test = "Chisq")$`Pr(>Chi)`[2],
      anova(unified_model_main, unified_model_score_by_distance, test = "Chisq")$`Pr(>Chi)`[2],
      anova(unified_model_main, unified_model_pam5_by_pam, test = "Chisq")$`Pr(>Chi)`[2]
    )
  ) %>%
  arrange(delta_aic, delta_bic)

selected_unified_model_name = interaction_model_summary$model[[1]]
selected_unified_model = unified_model_candidates[[selected_unified_model_name]]
selected_unified_model_table = tidy_unified_model(selected_unified_model) %>%
  mutate(selected_model = selected_unified_model_name)

selected_unified_model_focus_table = selected_unified_model_table %>%
  filter(feature_group %in% c("Guide score", "PAM upstream base", "PAM", "Edit position")) %>%
  mutate(
    term_order = c(
      "TTTC vs TTTA" = 1,
      "TTTG vs TTTA" = 2,
      "PAM 5' A vs C" = 3,
      "PAM 5' G vs C" = 4,
      "PAM 5' T vs C" = 5,
      "7-11 nt vs 1-6 nt" = 6,
      "12-17 nt vs 1-6 nt" = 7,
      "18-23 nt vs 1-6 nt" = 8,
      "PAM-overlap vs 1-6 nt" = 9,
      "DeepCpf1 score" = 10
    )[term_label]
  ) %>%
  arrange(term_order) %>%
  select(-term_order)

unified_model_data = unified_model_data %>%
  mutate(unified_model_prediction = predict(selected_unified_model, newdata = unified_model_data, type = "response"))

design_outcomes = design_outcomes %>%
  left_join(unified_model_data %>% select(uid, unified_model_prediction), by = "uid") %>%
  mutate(
    optimal_top_deepcpf1 = deepcpf1_quartile == "Q4 highest",
    optimal_pam_distance = !pam_dist_bin %in% c("PAM", "18-23"),
    optimal_non_tttg_pam = pam != "TTTG",
    optimal_non_c_pam5 = pam_upstream_nt != "C",
    optimal_all_filters = optimal_top_deepcpf1 & optimal_pam_distance & optimal_non_tttg_pam & optimal_non_c_pam5,
    optimal_filter_count = as.integer(optimal_top_deepcpf1) + as.integer(optimal_pam_distance) +
      as.integer(optimal_non_tttg_pam) + as.integer(optimal_non_c_pam5)
  )

optimal_filter_labels = c(
  optimal_top_deepcpf1 = "DeepCpf1 top quartile",
  optimal_pam_distance = "PAM distance not PAM/18-23",
  optimal_non_tttg_pam = "PAM not TTTG",
  optimal_non_c_pam5 = "5' PAM base not C"
)

optimal_filter_summary = bind_rows(
  summarize_design_subset(design_outcomes %>% filter(assayed_design), "All assayed designs"),
  summarize_design_subset(design_outcomes %>% filter(assayed_design, optimal_top_deepcpf1), "DeepCpf1 top quartile"),
  summarize_design_subset(design_outcomes %>% filter(assayed_design, optimal_pam_distance), "PAM distance not PAM/18-23"),
  summarize_design_subset(design_outcomes %>% filter(assayed_design, optimal_non_tttg_pam), "PAM not TTTG"),
  summarize_design_subset(design_outcomes %>% filter(assayed_design, optimal_non_c_pam5), "5' PAM base not C"),
  summarize_design_subset(design_outcomes %>% filter(assayed_design, optimal_all_filters), "All 4 optimal filters")
) %>%
  mutate(
    baseline_rate = design_success_rate[filter_label == "All assayed designs"][1],
    absolute_gain_vs_all = design_success_rate - baseline_rate
  )

three_way_filter_combinations = combn(names(optimal_filter_labels), 3, simplify = FALSE)

three_way_filter_summary = map_dfr(
  three_way_filter_combinations,
  function(filter_set) {
    filter_label = paste(unname(optimal_filter_labels[filter_set]), collapse = " + ")

    filtered_data = design_outcomes %>% filter(assayed_design)
    for (filter_name in filter_set) {
      filtered_data = filtered_data %>% filter(.data[[filter_name]])
    }

    summarize_design_subset(filtered_data, filter_label) %>%
      mutate(
        filter_group = "3-way combination",
        filter_set = paste(filter_set, collapse = "|")
      )
  }
) %>%
  mutate(
    baseline_rate = optimal_filter_summary$design_success_rate[optimal_filter_summary$filter_label == "All assayed designs"][1],
    absolute_gain_vs_all = design_success_rate - baseline_rate
  ) %>%
  arrange(desc(design_success_rate), desc(total_guides))

three_way_filter_plot_summary = bind_rows(
  optimal_filter_summary %>%
    filter(filter_label == "All assayed designs") %>%
    mutate(filter_group = "Baseline", filter_set = "baseline"),
  three_way_filter_summary,
  optimal_filter_summary %>%
    filter(filter_label == "All 4 optimal filters") %>%
    mutate(filter_group = "All 4 filters", filter_set = "all_4")
)

assayed_design_outcomes = design_outcomes %>%
  filter(assayed_design) %>%
  mutate(three_rule_set = optimal_pam_distance & optimal_non_tttg_pam & optimal_non_c_pam5)

set.seed(20250326)

three_rule_train_uids = assayed_design_outcomes %>%
  split(.$edited_design) %>%
  map(
    ~ sample(
      .x$uid,
      size = round(nrow(.x) * 0.7),
      replace = FALSE
    )
  ) %>%
  unlist(use.names = FALSE)

three_rule_train_test = assayed_design_outcomes %>%
  mutate(split = ifelse(uid %in% three_rule_train_uids, "Train", "Test"))

three_rule_train_test_summary = map_dfr(
  c("Train", "Test"),
  function(split_name) {
    split_data = three_rule_train_test %>% filter(split == split_name)

    bind_rows(
      summarize_design_subset(split_data, "All assayed designs"),
      summarize_design_subset(split_data %>% filter(three_rule_set), "3-rule set")
    ) %>%
      dplyr::rename(subset_label = filter_label) %>%
      mutate(split = split_name)
  }
) %>%
  group_by(split) %>%
  mutate(
    baseline_rate = design_success_rate[subset_label == "All assayed designs"][1],
    absolute_gain_vs_split_baseline = design_success_rate - baseline_rate,
    selected_fraction_of_split = total_guides / total_guides[subset_label == "All assayed designs"][1]
  ) %>%
  ungroup()

three_rule_train_test_metrics = three_rule_train_test %>%
  group_by(split) %>%
  summarize(
    total_guides = n(),
    selected_guides = sum(three_rule_set),
    selected_fraction = selected_guides / total_guides,
    true_positive = sum(three_rule_set & edited_design),
    false_positive = sum(three_rule_set & !edited_design),
    true_negative = sum(!three_rule_set & !edited_design),
    false_negative = sum(!three_rule_set & edited_design),
    precision = true_positive / (true_positive + false_positive),
    recall = true_positive / (true_positive + false_negative),
    specificity = true_negative / (true_negative + false_positive),
    negative_predictive_value = true_negative / (true_negative + false_negative),
    accuracy = (true_positive + true_negative) / total_guides,
    baseline_success_rate = mean(edited_design),
    absolute_gain_vs_baseline = precision - baseline_success_rate,
    enrichment_ratio_vs_baseline = precision / baseline_success_rate,
    .groups = "drop"
  )

selected_filter_comparison_summary = bind_rows(
  summarize_design_subset(design_outcomes %>% filter(assayed_design), "All assayed designs"),
  summarize_design_subset(
    design_outcomes %>% filter(assayed_design, optimal_pam_distance, optimal_non_tttg_pam, optimal_non_c_pam5),
    "PAM distance not PAM/18-23 + PAM not TTTG + 5' PAM base not C"
  ),
  summarize_design_subset(
    design_outcomes %>% filter(assayed_design, optimal_top_deepcpf1, optimal_pam_distance, optimal_non_tttg_pam),
    "DeepCpf1 top quartile + PAM distance not PAM/18-23 + PAM not TTTG"
  ),
  summarize_design_subset(
    design_outcomes %>% filter(assayed_design, optimal_all_filters),
    "All 4 optimal filters"
  )
) %>%
  mutate(
    baseline_rate = design_success_rate[filter_label == "All assayed designs"][1],
    absolute_gain_vs_all = design_success_rate - baseline_rate
  )

designs_assayed_optimal = design_outcomes %>%
  filter(
    assayed_design,
    optimal_all_filters
  ) %>%
  arrange(desc(unified_model_prediction), uid)

## plots from the notebook
pam_distance_heatmap = pam_heatmap_summary %>%
  mutate(
    pam = factor(pam, levels = c("TTTC", "TTTA", "TTTG")),
    pam_dist_bin = factor(pam_dist_bin, levels = c("PAM", "1-6", "7-11", "12-17", "18-23")),
    label = sprintf("%s\nn=%d", percent(design_success_rate, accuracy = 1), total_guides)
  ) %>%
  ggplot(aes(x = pam_dist_bin, y = pam, fill = design_success_rate)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradientn(colors = c("white", "#E69F00"), labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Edit position relative to PAM", y = "PAM", fill = "Design success") +
  theme_design_rules(aspect_ratio = 0.8, legend_position = "right")

deepcpf1_plot = plot_rate_bars(deepcpf1_summary, "deepcpf1_quartile", "DeepCpf1 quartile", "Higher DeepCpf1 scores enrich successful designs")
pam_upstream_plot = plot_rate_bars(pam_upstream_summary, "pam_upstream_nt", "5' PAM base (N in NTTTV)", "The upstream PAM base adds additional design signal")
rule_count_plot = plot_rate_bars(rule_count_summary, "preferred_rule_count", "Preferred-rule count", "Simple rule bundles separate strong from weak designs")

multi_colony_consistency_plot_data = multi_colony_consistency_summary %>%
  mutate(
    assayed_colonies = factor(assayed_colonies),
    label = paste0("n=", total_designs)
  )

multi_colony_consistency_label_data = multi_colony_consistency_plot_data %>%
  distinct(assayed_colonies, total_designs, label)

multi_colony_consistency_plot = ggplot(multi_colony_consistency_plot_data, aes(x = assayed_colonies, y = fraction, fill = consistency)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  geom_text(
    data = multi_colony_consistency_label_data,
    aes(x = assayed_colonies, y = 1.03, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_fill_manual(values = c("All edited" = "#355070", "Mixed" = "#E56B6F", "None edited" = "#D9D9D9")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.08), expand = c(0, 0)) +
  labs(x = "Assayed colonies per design", y = "Fraction of designs", fill = NULL, title = "Repeated colonies do not always agree") +
  theme_design_rules(aspect_ratio = 1, legend_position = "top")

dr_perturbed_plot = plot_rate_bars(dr_perturbed_summary, "dr_perturbed_flag", "DR folding call", "Perturbed DR folding alone was not a major failure mode")
dr_canonical_plot = plot_rate_bars(dr_canonical_summary, "dr_canonical_flag", "DR structure class", "Canonical DR structures were not strongly enriched")
dr_spacer_pair_plot = plot_rate_bars(dr_spacer_pair_summary, "dr_spacer_pair_flag", "DR-spacer pairing", "Any DR-spacer pairing had only a weak marginal effect")

guide_gc_23_plot = plot_rate_bars(guide_gc_23_summary, "guide_gc_23_bin", "Guide GC (23 nt)", "Full spacer GC content")
guide_gc_20_plot = plot_rate_bars(guide_gc_20_summary, "guide_gc_20_bin", "Guide GC (first 20 nt)", "Guide GC, first 20 nt")
guide_gc_18_plot = plot_rate_bars(guide_gc_18_summary, "guide_gc_18_bin", "Guide GC (first 18 nt)", "Guide GC, first 18 nt")
donor_wt_gc_plot = plot_rate_bars(donor_wt_gc_summary, "donor_wt_gc_bin", "WT donor GC", "WT donor GC content")

attrition_plot = stage_summary %>%
  filter(round != "overall") %>%
  transmute(
    round = round,
    picked_designs = picked_designs,
    assayed_designs = assayed_designs,
    edited_designs = edited_designs,
    picked_colonies = picked_colonies,
    assayed_colonies = assayed_colonies,
    edited_colonies = edited_colonies
  ) %>%
  pivot_longer(cols = -round, names_to = "metric", values_to = "count") %>%
  mutate(
    measure = ifelse(str_detect(metric, "design"), "Designs", "Colonies"),
    stage = case_when(
      str_starts(metric, "picked") ~ "Picked",
      str_starts(metric, "assayed") ~ "Assayed",
      TRUE ~ "Edited"
    ),
    stage = factor(stage, levels = c("Picked", "Assayed", "Edited")),
    round = recode(round, round_1 = "Round 1", round_2 = "Round 2")
  ) %>%
  ggplot(aes(x = stage, y = count, fill = round)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  facet_wrap(~measure, scales = "free_y") +
  scale_fill_manual(values = c("Round 1" = "#6D597A", "Round 2" = "#B56576")) +
  labs(x = NULL, y = "Count", fill = NULL, title = "Genome-wide assay attrition") +
  theme_design_rules(aspect_ratio = NULL, legend_position = "top")

coverage_qc_plot = ggplot(designs_all_annotated %>% filter(!is.na(mean_cov), !is.na(target_basecov)), aes(x = mean_cov, y = target_basecov, color = as.factor(edited_colony == 1))) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey50") +
  geom_smooth(method = "lm", se = FALSE, color = "#355070", linewidth = 0.8) +
  scale_color_manual(values = c("TRUE" = "#B56576", "FALSE" = "#6C757D"), labels = c("FALSE" = "Not edited", "TRUE" = "Edited")) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.2,
    label = sprintf("Pearson r = %.2f", suppressWarnings(cor(designs_all_annotated$mean_cov, designs_all_annotated$target_basecov, use = "pairwise.complete.obs"))),
    size = 3.4
  ) +
  labs(x = "Mean genome coverage", y = "Target base coverage", color = NULL, title = "Coverage QC") +
  theme_design_rules(aspect_ratio = 1, legend_position = "top")

max_A_plot = plot_rate_bars(max_A_summary, "max_A", "max A run", "Long A-runs are at most a secondary filter")
max_T_plot = plot_rate_bars(max_T_summary, "max_T", "max T run", "T-runs were not a clean marginal failure mode")
region_class_plot = plot_rate_bars(region_class_summary, "region_class", "Variant region class", "Locus context affects editability but is not a cassette rule")
hotspot_plot = plot_rate_bars(hotspot_summary, "hotspot_flag", "SV hotspot annotation", "SV hotspots are a contextual caution flag")

unified_model_plot = plot_model_odds_ratios(
  selected_unified_model_table,
  paste0("Unified binary cassette-level model: ", selected_unified_model_name)
)
unified_model_focus_plot = plot_model_odds_ratios(
  selected_unified_model_focus_table,
  paste0("Cassette editing odds for key design parameters: ", selected_unified_model_name)
)

optimal_filter_plot = plot_optimal_filter_effects(optimal_filter_summary)
three_way_filter_plot = plot_three_way_filter_effects(three_way_filter_plot_summary)
three_rule_train_test_plot = plot_rule_train_test_effects(three_rule_train_test_summary)
selected_filter_comparison_plot = plot_selected_filter_comparison(selected_filter_comparison_summary)
chromosome_circos_plot = plot_chromosome_editing_circos(chromosome_circos_summary, chromosome_layout, chromosome_bin_width)
stacked_design_outcome_core_plot = plot_stacked_design_outcomes(
  stacked_design_outcome_summary %>% filter(parameter_group == "Core cassette parameters"),
  "Design outcomes across core cassette parameters",
  "Designs with any high-confidence non-target variant are grouped under 'Other variants observed'",
  ncol = 2
)
stacked_design_outcome_extended_plot = plot_stacked_design_outcomes(
  stacked_design_outcome_summary %>% filter(parameter_group == "GC and DR parameters"),
  "Design outcomes across GC and DR-related parameters",
  "Assayed designs only; bars show intended edit only, other variants observed, or unedited",
  ncol = 2
)

background_variant_burden_plot = background_colony_summary %>%
  mutate(target_edit_status = ifelse(edited_colony, "Edited target", "Unedited target")) %>%
  ggplot(aes(x = target_edit_status, y = high_conf_non_target_count, fill = target_edit_status)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.14, height = 0.02, alpha = 0.45, size = 1.1) +
  scale_fill_manual(values = c("Edited target" = "#355070", "Unedited target" = "#B56576")) +
  scale_y_continuous(breaks = 0:4, limits = c(0, 4)) +
  labs(x = NULL, y = "High-confidence non-target variants per colony", fill = NULL, title = "Background variant burden is sparse") +
  theme_design_rules(aspect_ratio = 1, legend_position = "none")

candidate_main_figure = ((pam_distance_heatmap | deepcpf1_plot) / (pam_upstream_plot | multi_colony_consistency_plot)) +
  plot_annotation(tag_levels = "A")

dr_support_figure = (pam_upstream_plot | dr_perturbed_plot | dr_spacer_pair_plot) +
  plot_annotation(tag_levels = "A")

gc_dependency_figure = ((guide_gc_23_plot | guide_gc_20_plot) / (guide_gc_18_plot | donor_wt_gc_plot)) +
  plot_annotation(tag_levels = "A")

supplementary_figure = ((attrition_plot | coverage_qc_plot) / (max_A_plot | max_T_plot) / (region_class_plot | hotspot_plot)) +
  plot_annotation(tag_levels = "A")

## save plots
save_plot(pam_distance_heatmap, file.path(plots_dir, "pam_distance_heatmap.svg"), width = 6.5, height = 4.5)
save_plot(deepcpf1_plot, file.path(plots_dir, "deepcpf1_quartiles.pdf"), width = 5.2, height = 4.8)
save_plot(pam_upstream_plot, file.path(plots_dir, "pam_upstream_nt.pdf"), width = 5.2, height = 4.8)
save_plot(rule_count_plot, file.path(plots_dir, "preferred_rule_count.pdf"), width = 5.2, height = 4.8)
save_plot(multi_colony_consistency_plot, file.path(plots_dir, "multi_colony_consistency.pdf"), width = 5.5, height = 5)
save_plot(dr_perturbed_plot, file.path(plots_dir, "dr_perturbed_summary.pdf"), width = 5.4, height = 4.8)
save_plot(dr_canonical_plot, file.path(plots_dir, "dr_canonical_summary.pdf"), width = 5.4, height = 4.8)
save_plot(dr_spacer_pair_plot, file.path(plots_dir, "dr_spacer_pair_summary.pdf"), width = 5.4, height = 4.8)
save_plot(guide_gc_23_plot, file.path(plots_dir, "guide_gc_23_summary.pdf"), width = 5.2, height = 4.8)
save_plot(guide_gc_20_plot, file.path(plots_dir, "guide_gc_20_summary.pdf"), width = 5.2, height = 4.8)
save_plot(guide_gc_18_plot, file.path(plots_dir, "guide_gc_18_summary.pdf"), width = 5.2, height = 4.8)
save_plot(donor_wt_gc_plot, file.path(plots_dir, "donor_wt_gc_summary.pdf"), width = 5.2, height = 4.8)
save_plot(attrition_plot, file.path(plots_dir, "attrition_plot.pdf"), width = 8, height = 4.5)
save_plot(coverage_qc_plot, file.path(plots_dir, "coverage_qc.pdf"), width = 5.5, height = 5.5)
save_plot(max_A_plot, file.path(plots_dir, "max_A_summary.pdf"), width = 5.2, height = 4.8)
save_plot(max_T_plot, file.path(plots_dir, "max_T_summary.pdf"), width = 5.2, height = 4.8)
save_plot(region_class_plot, file.path(plots_dir, "region_class_summary.pdf"), width = 5.6, height = 4.8)
save_plot(hotspot_plot, file.path(plots_dir, "hotspot_summary.pdf"), width = 5.2, height = 4.8)
save_plot(unified_model_plot, file.path(plots_dir, "unified_design_model_odds_ratios.pdf"), width = 9, height = 5.5)
save_plot(unified_model_focus_plot, file.path(plots_dir, "unified_design_model_key_parameters.pdf"), width = 8.8, height = 5.4)
save_plot(optimal_filter_plot, file.path(plots_dir, "optimal_filter_effects.pdf"), width = 8.8, height = 5.8)
save_plot(three_way_filter_plot, file.path(plots_dir, "three_way_filter_effects.pdf"), width = 9.8, height = 6.2)
save_plot(three_rule_train_test_plot, file.path(plots_dir, "three_rule_train_test_validation.pdf"), width = 9.6, height = 5.6)
save_plot(selected_filter_comparison_plot, file.path(plots_dir, "selected_filter_comparison_orange.svg"), width = 9.2, height = 5.8)
save_plot(chromosome_circos_plot, file.path(plots_dir, "chromosome_editing_circos.pdf"), width = 9.6, height = 9.6)
save_plot(stacked_design_outcome_core_plot, file.path(plots_dir, "stacked_design_outcomes_core_parameters.pdf"), width = 12.5, height = 9.5)
save_plot(stacked_design_outcome_extended_plot, file.path(plots_dir, "stacked_design_outcomes_extended_parameters.pdf"), width = 12.5, height = 11)
save_plot(background_variant_burden_plot, file.path(plots_dir, "background_variant_burden.pdf"), width = 5.5, height = 5)
save_plot(candidate_main_figure, file.path(plots_dir, "candidate_main_figure_design_rules.pdf"), width = 11, height = 10)
save_plot(dr_support_figure, file.path(plots_dir, "pam_upstream_and_dr_folding.pdf"), width = 15, height = 4.8)
save_plot(gc_dependency_figure, file.path(plots_dir, "gc_dependency_panels.pdf"), width = 11, height = 9.5)
save_plot(supplementary_figure, file.path(plots_dir, "supplementary_design_rule_context.pdf"), width = 11, height = 13)

## save tables
write.csv(designs_all_annotated, file.path(tables_dir, "colony_assays.csv"), row.names = FALSE)
write.csv(design_outcomes, file.path(tables_dir, "design_outcomes.csv"), row.names = FALSE)
write.csv(stage_summary, file.path(tables_dir, "stage_summary.csv"), row.names = FALSE)
write.csv(pam_summary, file.path(tables_dir, "pam_summary.csv"), row.names = FALSE)
write.csv(pam_distance_summary, file.path(tables_dir, "pam_distance_summary.csv"), row.names = FALSE)
write.csv(pam_heatmap_summary, file.path(tables_dir, "pam_distance_heatmap_summary.csv"), row.names = FALSE)
write.csv(deepcpf1_summary, file.path(tables_dir, "deepcpf1_summary.csv"), row.names = FALSE)
write.csv(pam_upstream_summary, file.path(tables_dir, "pam_upstream_summary.csv"), row.names = FALSE)
write.csv(full_pam_5mer_summary, file.path(tables_dir, "full_pam_5mer_summary.csv"), row.names = FALSE)
write.csv(ntttv_pattern_summary, file.path(tables_dir, "ntttv_pattern_summary.csv"), row.names = FALSE)
write.csv(rule_count_summary, file.path(tables_dir, "preferred_rule_count_summary.csv"), row.names = FALSE)
write.csv(dr_perturbed_summary, file.path(tables_dir, "dr_perturbed_summary.csv"), row.names = FALSE)
write.csv(dr_canonical_summary, file.path(tables_dir, "dr_canonical_summary.csv"), row.names = FALSE)
write.csv(dr_spacer_pair_summary, file.path(tables_dir, "dr_spacer_pair_summary.csv"), row.names = FALSE)
write.csv(guide_gc_23_summary, file.path(tables_dir, "guide_gc_23_summary.csv"), row.names = FALSE)
write.csv(guide_gc_20_summary, file.path(tables_dir, "guide_gc_20_summary.csv"), row.names = FALSE)
write.csv(guide_gc_18_summary, file.path(tables_dir, "guide_gc_18_summary.csv"), row.names = FALSE)
write.csv(donor_wt_gc_summary, file.path(tables_dir, "donor_wt_gc_summary.csv"), row.names = FALSE)
write.csv(max_A_summary, file.path(tables_dir, "max_A_summary.csv"), row.names = FALSE)
write.csv(max_T_summary, file.path(tables_dir, "max_T_summary.csv"), row.names = FALSE)
write.csv(region_class_summary, file.path(tables_dir, "region_class_summary.csv"), row.names = FALSE)
write.csv(hotspot_summary, file.path(tables_dir, "hotspot_summary.csv"), row.names = FALSE)
write.csv(multi_colony_consistency_summary, file.path(tables_dir, "multi_colony_consistency_summary.csv"), row.names = FALSE)
write.csv(design_rule_summary, file.path(tables_dir, "design_rule_summary.csv"), row.names = FALSE)
write.csv(background_variants, file.path(tables_dir, "all_variants_with_score_sv.csv"), row.names = FALSE)
write.csv(high_conf_non_target_variants, file.path(tables_dir, "high_conf_non_target_variants.csv"), row.names = FALSE)
write.csv(background_colony_summary, file.path(tables_dir, "background_colony_summary.csv"), row.names = FALSE)
write.csv(background_edit_status_summary, file.path(tables_dir, "background_edit_status_summary.csv"), row.names = FALSE)
write.csv(recurrent_background_variants, file.path(tables_dir, "recurrent_background_variants.csv"), row.names = FALSE)
write.csv(interaction_model_summary, file.path(tables_dir, "interaction_model_comparison.csv"), row.names = FALSE)
write.csv(selected_unified_model_table, file.path(tables_dir, "unified_design_model_coefficients.csv"), row.names = FALSE)
write.csv(selected_unified_model_focus_table, file.path(tables_dir, "unified_design_model_key_parameters.csv"), row.names = FALSE)
write.csv(chromosome_circos_summary, file.path(tables_dir, "chromosome_editing_circos_summary.csv"), row.names = FALSE)
write.csv(stacked_design_outcome_summary, file.path(tables_dir, "stacked_design_outcome_summary.csv"), row.names = FALSE)
write.csv(optimal_filter_summary, file.path(tables_dir, "optimal_filter_summary.csv"), row.names = FALSE)
write.csv(three_way_filter_summary, file.path(tables_dir, "three_way_filter_summary.csv"), row.names = FALSE)
write.csv(three_rule_train_test_summary, file.path(tables_dir, "three_rule_train_test_summary.csv"), row.names = FALSE)
write.csv(three_rule_train_test_metrics, file.path(tables_dir, "three_rule_train_test_metrics.csv"), row.names = FALSE)
write.csv(selected_filter_comparison_summary, file.path(tables_dir, "selected_filter_comparison_summary.csv"), row.names = FALSE)
write.csv(designs_assayed_optimal, file.path(tables_dir, "designs_assayed_optimal.csv"), row.names = FALSE)
write.csv(deepcpf1_prediction_summary, file.path(tables_dir, "deepcpf1_prediction_summary.csv"), row.names = FALSE)

print(tibble(
  metric = c(
    "Picked colonies",
    "Assayed colonies",
    "Edited colonies",
    "Unique picked designs",
    "Assayed designs",
    "Successful designs",
    "Direct-repeat sequence",
    "High-confidence non-target variants"
  ),
  value = c(
    nrow(designs_all_annotated),
    sum(designs_all_annotated$assayed_colony),
    sum(designs_all_annotated$edited_colony),
    nrow(design_outcomes),
    sum(design_outcomes$assayed_design),
    sum(design_outcomes$edited_design),
    DR_seq,
    nrow(high_conf_non_target_variants)
  )
))


## pull out designs that have a mean coverage of at least 5 that are in hotspots
