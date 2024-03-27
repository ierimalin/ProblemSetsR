# Code to generate tables of means and compares them for tests

library(tidyverse)
library(AER)
library(kableExtra)
library(stargazer)
library(lubridate)
library(stringr)

options(knitr.kable.NA = '')

                                        # Load data

df_all <- read_csv("clean_data.csv") %>%
  mutate(
    gender    = as_factor(gender),
    treatment = fct_relevel(as_factor(treatment), "BROAD", "NARROW", "LOW", "PARTIAL", "BEFORE", "AFTER"),
    scenario  = as_factor(scenario)
  ) %>%
  ## Rename the treatments in line with the text
  mutate(
    treatment = dplyr::recode(treatment,
                              BROAD   = "NONE",
                              NARROW  = "BOTH",
                              LOW     = "MONEY/LOW",
                              PARTIAL = "MONEY")
  )

consistent_df_all <- df_all %>%
  ## Keep only choices which were done consistently (drop only that scenario)
  filter(inconsistent_choices == 0) %>%
  ## Keep only people who finish the study and hence the tasks
  filter(!is.na(time_questionnaire2)) %>%
  select(pid, scenario, treatment, reservation_wage, inconsistent_choices, everything())

## To avoid reusing inconsistent data by mistake
df_all <- NULL

consistent_dfba <- consistent_df_all %>%
  filter(treatment %in% c("BOTH", "BEFORE", "AFTER", "NONE"))
consistent_df <- consistent_df_all %>%
  filter(treatment %in% c("MONEY/LOW", "BOTH", "MONEY", "NONE"))

# Bar and kernel density plots
raw_bar_plot <- ggplot(consistent_df, aes(x = reservation_wage)) +
  geom_bar(aes(y = ..prop..)) +
  facet_grid(treatment ~ scenario) +
  labs(x = "Reservation Wage", y = "Proportion")

ggsave("bar_plot.png", plot = raw_bar_plot, "png")

kernel_density_plot <- ggplot(consistent_df, aes(x = reservation_wage, group = treatment, color = treatment)) +
  geom_density() +
  facet_wrap(~ scenario) +
  labs(x = "Reservation Wage", y = "Density")

ggsave("density_plot.png", plot = kernel_density_plot, "png")

# Various outputs and tests about means
means_data <- function(df) {
  paste_to_percent <- function(v) paste0(as.character(round(v)), "%")
  df %>%
    filter( !is.na(reservation_wage)) %>%
    group_by(scenario, treatment) %>%
    summarise(
      res_wage   = mean(reservation_wage),
      standard_error = sd(reservation_wage)/sqrt(n()),
      share_upper_bound = paste_to_percent(100 * mean(reservation_wage > 4.01)),
      N = n()
    ) %>%
    ungroup()
}

table_of_means_data <- function(means_df, caption, label, ntreatments = 4) {
  kbl(
    means_df,
    "latex",
    col.names = c("Treatment", "Res. Wage", "Std Err", "% upper bound", "N"),
    booktabs = T,
    caption = caption,
    align = 'lrrrr',
    digits=2,
    label=label
  ) %>%
    kable_styling(latex_options=c("hold_position"), font_size=12) %>%
    pack_rows("Scenario 1", 1, ntreatments) %>%
    pack_rows("Scenario 2", ntreatments + 1, 2*ntreatments)
}

main_treatments <- c("NONE", "BOTH", "MONEY/LOW")
main_df <- consistent_df %>%
  filter(treatment %in% main_treatments)

partial_treatments <- c("NONE", "BOTH", "MONEY")
partial_df <- consistent_df %>%
  filter(treatment %in% partial_treatments)

main_df %>%
  means_data() %>%
  select(-scenario) %>%
  table_of_means_data(
    "Means of main treatments by scenario",
    label="means_data",
    ntreatments=length(main_treatments)) %>%
  write("means_data.tex")

january1 <- ymd("20200101")

consistent_df_presentation_on_page <- consistent_df %>%
  filter(!(treatment == "BOTH" & date(session_date) <= january1))

consistent_df_presentation_on_page %>%
  means_data() %>%
  select(-scenario) %>%
  table_of_means_data(
    caption="Means of main treatments by scenario, restricting BOTH to sessions where baseline is revealed on first choice page.",
    label="means_data_on_page"
  ) %>%
  write("means_presentation_on_page.tex")

# Same information when limiting to data with treatment BOTH having endowment
# revealed right before the first choice page
consistent_df_presentation_before_page <- consistent_df %>%
  filter(!(treatment == "BOTH" & date(session_date) >= january1))

consistent_df_presentation_before_page %>%
  means_data() %>%
  select(-scenario) %>%
  table_of_means_data(
    caption="Means of main treatments by scenario, restricting BOTH to sessions where baseline is revealed right before the first choice page.",
    label="means_data_before_page"
  ) %>%
  write("means_presentation_before_page.tex")

# t-tests with clustered standard errors (aka OLS with clustering) when
# pooling scenarios
lm_pooled <- lm(reservation_wage ~ scenario + treatment, data = consistent_df)
lm_by_scenario <- lm(reservation_wage ~ scenario*treatment, data = consistent_df)

get_clustered_se <- function(m, df) {
  vcov_cluster <- multiwayvcov::cluster.vcov(m, df$pid)
  coef_clustered <- coeftest(m, vcov_cluster)
}

clustered <- get_clustered_se(lm_pooled, consistent_df)

stargazer(lm_pooled, lm_pooled,
          title="Linear Regressions of reservation wages by treatment averaged across scenarios. With and without clustered standard errors by participant. The default treatment is NONE, so that broad bracketing predicts a null estimate for the fixed effect for BOTH (rejected), and narrow bracketing predicts equal fixed effects for BOTH and MONEY/FALSE (not rejected). Broad bracketing of money (but not necessarily work) predicts a null estimate for MONEY (rejected); broad bracketing of work (but not necessarily money) predicts equal fixed effects for BOTH and MONEY (not rejected). The last-mentioned test is not rejected because the difference is positive in Scenario 1 and negative in Scenario 2, which is average in this regression. Since bracketing makes predictions that should hold for every scenario, the test in our main text is the better test (accounting for double-testing).",
          header=FALSE,
          omit.stat=c("adj.rsq", "ser", "ll", "wald"),
          dep.var.labels="Reservation wage",
          column.labels=c("No clustering", "Clustering by participant"),
          model.names=FALSE,
          star.char=c("", "", ""),
          notes = c("Standard errors in parentheses."),
          digits = 2,
          label = "tab:linear-regressions",
          notes.append = FALSE,
          se = list(NULL, clustered[,2]),
          t = list(NULL, NULL, NULL, clustered[,3]),
          p = list(NULL, NULL, NULL, clustered[,4])
          ) %>%
  write("linear_regression_pooled_and_by_scenario.tex")

## Tables of p-values

pvalue_to_char <- function(p) {
  ifelse(p < 0.001, "$< 0.001$", paste0("$", sprintf(p, fmt = '%#.3f'), "$"))
}

wilcox_pvalue <- function(x,y) {
  oldWarn = getOption("warn")
  options(warn=-1)
  res <- pvalue_to_char(wilcox.test(x, y, exact=TRUE)$p.value)
  options(warn=oldWarn)
  res
}

t_pvalue <- function(x,y) pvalue_to_char(t.test(x,y)$p.value)

get_hyp_pvalue <- function(df, f_pvalue, treatments) {
  l <- list()
  for (t in treatments) {
    l[[t]] <- df$reservation_wage[df$treatment == t]
  }

  n <- length(treatments)
  pvalues <- list()
  for (i in 1:(n - 1)) {
    pvalues[[ treatments[i] ]] <- c()
    for (j in 2:n) {
      if (j > i) {
        pvalues[[ treatments[i] ]] <- c(
          pvalues[[ treatments[i] ]],
          f_pvalue(l[[ treatments[i] ]], l[[ treatments[j] ]])
        )
      } else {
        pvalues[[ treatments[i] ]] <- c(pvalues[[ treatments[i] ]], NA)
      }
    }
  }

  result_table <- tibble(Treatments = treatments[2:n])

  for ( i in 1:(n - 1) ) {
    result_table[[ treatments[i] ]] <- pvalues[[ treatments[i] ]]
  }
  result_table
}

get_scenario_tables <- function(df, ptest, treatments) {
  table1_presentation <- get_hyp_pvalue(filter(df, scenario == "Scenario1"), ptest, treatments)
  table2_presentation <- get_hyp_pvalue(filter(df, scenario == "Scenario2"), ptest, treatments)
  rbind( table1_presentation, table2_presentation )
}

tex_for_pvalue <- function(df, ptest, filename, caption, label, treatments) {
  t <- get_scenario_tables(df, ptest, treatments)
  n <- length(treatments)
  kbl(
    t,
    "latex",
    escape=FALSE,
    booktabs=T,
    caption=caption,
    align="lccc",
    label=label,
  ) %>%
    kable_styling(font_size = 12) %>%
    pack_rows("Scenario 1", 1, n - 1) %>%
    pack_rows("Scenario 2", n, 2*(n - 1)) %>%
    write(filename)
}

tex_for_pvalue_all_plus_gender <- function(df, ptest, filename, caption, label, treatments) {
  t_all <- get_scenario_tables(df, ptest, treatments)
  t_female <- get_scenario_tables(df %>% filter( gender == "Female"), ptest, treatments) %>%
    select(-Treatments)
  t_male <- get_scenario_tables(df %>% filter( gender == "Male"), ptest, treatments) %>%
    select(-Treatments)

  for (coln in colnames(t_female)) {
    t_female[[ paste0(coln, "/F") ]] <- t_female[[ coln ]]
    t_female[[ coln ]] <- NULL
  }

  for (coln in colnames(t_male)) {
    t_male[[ paste0(coln, "/M") ]] <- t_male[[ coln ]]
    t_male[[ coln ]] <- NULL
  }

  t <- cbind(t_all, t_female, t_male)
  n <- length(treatments)
  kbl(
    t,
    "latex",
    escape=FALSE,
    booktabs=T,
    caption=caption,
    align="lcccccc",
    label=label,
  ) %>%
    kable_styling(font_size = 12) %>%
    pack_rows("Scenario 1", 1, n - 1) %>%
    pack_rows("Scenario 2", n, 2*(n - 1)) %>%
    add_header_above(c("Pooled" = 3, "Female" = 2, "Male" = 2)) %>%
    write(filename)
}

tex_for_pvalue_all_plus_gender(
  main_df,
  wilcox_pvalue,
  caption="Between-treatment p-values for main treatments based on two-sided Wilcoxon rank-sum tests, treating each individual in each scenario as a single independent observation. The first two columns are for pooled data, the next two for data restricted to female participants, the final two for data restricted to male participants.",
  filename="mwu_all_tables.tex",
  label="mwu_all_tables",
  treatments=main_treatments
)

tex_for_pvalue_all_plus_gender(
  main_df,
  t_pvalue,
  caption="Between-treatment p-values for main treatments based on two-sided t-tests by scenario, treating each individual in each scenario as a single independent observation. The first two columns are for pooled data, the next two for data restricted to female participants, the final two for data restricted to male participants.",
  filename="t_test_all_tables.tex",
  label="t_test_all_tables",
  treatments=main_treatments
)

tex_for_pvalue(
  main_df,
  wilcox_pvalue,
  caption="Between-treatment p-values for main treatments based on two-sided Wilcoxon rank-sum tests, treating each individual in each scenario as a single independent observation.",
  filename="mwu_tables.tex",
  label="mwu_tables",
  treatments=main_treatments
)

tex_for_pvalue(
  consistent_df,
  t_pvalue,
  caption="Between-treatment p-values for main treatments based on two-sided t-tests by scenario",
  filename="t_test_tables.tex",
  label="t_tables",
  treatments=main_treatments
)

tex_for_pvalue(
  partial_df,
  t_pvalue,
  caption="Between-treatment p-values for treatments NONE, MONEY, and BOTH based on two-sided t-tests, treating each individual in each scenario as a single independent observation.",
  filename="t_test_tables_partial.tex",
  label="t_test_tables_partial",
  treatments=partial_treatments
)

## 'clustered t-tests' that is, OLS pooling the data between scenarios and
## clustering at the individual level

get_clustered_se <- function(m, df) {
  vcov_cluster <- multiwayvcov::cluster.vcov(m, df$pid)
  coef_clustered <- coeftest(m, vcov_cluster)
}

lm_cluster_pvalue <- function(df, t1, t2) {
  d <- filter(df, treatment %in% c(t1, t2))
  lm_pooled <- lm(reservation_wage ~ scenario + treatment, data = d)
  clustered <- get_clustered_se(lm_pooled, d)
  clustered[3,4] # FIXME: brittle to refer to the pvalue by row and column
}

cluster_table <- function(df) {
  f <- function(t1, t2) {
    pvalue_to_char(lm_cluster_pvalue(df, t1, t2))
  }
  pvalue_none <- c(f("NONE", "BOTH"), f("NONE", "MONEY/LOW"), f("NONE", "MONEY"))
  pvalue_both <- c(NA , f("BOTH", "MONEY/LOW"), f("BOTH", "MONEY"))
  pvalue_low <- c(NA, NA, f("MONEY/LOW", "MONEY"))
  (pvalue_table1 <- tibble(Treatments = c("BOTH", "MONEY/LOW", "MONEY"), NONE = pvalue_none, BOTH = pvalue_both, `MONEY/LOW` = pvalue_low))
}

cluster_table_before_after <- function(df) {
  f <- function(t1, t2) {
    pvalue_to_char(lm_cluster_pvalue(df, t1, t2))
  }
  pvalue_both <- c(f("BEFORE", "BOTH"), f("BOTH", "AFTER"))
  pvalue_before <- c(NA, f("BEFORE", "AFTER"))
  (pvalue_table1 <- tibble(Treatments = c("BEFORE", "AFTER"), BOTH = pvalue_both, BEFORE = pvalue_before))
}

kbl(
  cluster_table(consistent_df),
  "latex",
  escape=FALSE,
  booktabs=T,
  caption="Between-treatment p-values for differences in means between treatments, averaging across scenarios: differences of treatment fixed-effects clustered at the individual level.",
  align="lccc",
  label = "clustered_FE_table"
) %>%
  kable_styling(font_size = 12) %>%
  write("clustered_fixed_effects_table.tex")

kbl(
  cluster_table_before_after(consistent_dfba),
  "latex",
  escape=FALSE,
  booktabs=T,
  caption="Between-treatment p-values for differences in means between treatments BEFORE, AFTER, and BOTH; averaging across scenarios: differences of treatment fixed-effects clustered at the individual level.",
  align="lccc",
  label = "clustered_FE_table_before_after"
) %>%
  kable_styling(font_size = 12) %>%
  write("clustered_fixed_effects_table_before_after.tex")

## Comparing variations of BOTH with other treatments and itself, depending on
## whether BOTH displayed endowment on the page before the first choice page or
## only on the first choice page.

tex_for_pvalue(
  consistent_df_presentation_on_page,
  ptest=t_pvalue,
  filename="t_test_presentation_on_page.tex",
  caption="Between-treatment p-values for main treatments based on two-sided t-tests by scenario. Restricted to those sessions of BOTH where baseline is revealed only on first choice page.",
  label='t_test_on_page',
  treatments=main_treatments
)

tex_for_pvalue(
  consistent_df_presentation_before_page,
  ptest=t_pvalue,
  filename="t_test_presentation_before_page.tex",
  caption="Between-treatment p-values for main treatments based on two-sided t-tests by scenario. Restricted to those sessions of BOTH where baseline is revealed right before the first choice page.",
  label='t_test_before_page',
  treatments=main_treatments
)

tex_for_pvalue(
  consistent_df_presentation_on_page,
  ptest=wilcox_pvalue,
  filename="mwu_presentation_on_page.tex",
  caption="Between-treatment p-values for main treatments based on two-sided Wilcoxon rank-sum tests, treating each individual in each scenario as a single independent observation. Restricted to those sessions of BOTH where baseline is revealed only on first choice page.",
  label='mwu_test_on_page',
  treatments=main_treatments
)

tex_for_pvalue(
  consistent_df_presentation_on_page,
  ptest=t_pvalue,
  filename="t_test_presentation_on_page.tex",
  caption="Between-treatment p-values for main treatments based on two-sided t-test, treating each individual in each scenario as a single independent observation. Restricted to those sessions of BOTH where baseline is revealed only on first choice page.",
  label='t_test_on_page',
  treatments=main_treatments
)

tex_for_pvalue(
  consistent_df_presentation_before_page,
  ptest=wilcox_pvalue,
  filename="mwu_presentation_before_page.tex",
  caption="Between-treatment p-values for main treatments based on two-sided Wilcoxon rank-sum tests, treating each individual in each scenario as a single independent observation. Restricted to those sessions of BOTH where baseline is revealed right before the first choice page.",
  label='mwu_test_before_page',
  treatments=main_treatments
)

tex_for_pvalue(
  consistent_df_presentation_before_page,
  ptest=t_pvalue,
  filename="t_test_presentation_before_page.tex",
  caption="Between-treatment p-values for main treatments based on two-sided t-test, treating each individual in each scenario as a single independent observation. Restricted to those sessions of BOTH where baseline is revealed right before the first choice page.",
  label='t_test_test_before_page',
  treatments=main_treatments
)



## BOTH with display of endowment before and on first choice page

tex_for_narrow_v_narrow <- function(ptest, ptestlabel, caption, filename, label) {
  narrow_before_page <- consistent_df_presentation_before_page %>%
    filter(treatment == "BOTH")
  narrow_on_page <- consistent_df_presentation_on_page %>%
    filter(treatment == "BOTH")

  nvn1 <- ptest(
    narrow_before_page$reservation_wage[narrow_before_page$scenario == "Scenario1"],
    narrow_on_page$reservation_wage[narrow_before_page$scenario == "Scenario1"]
  )

  nvn2 <- ptest(
    narrow_before_page$reservation_wage[narrow_before_page$scenario == "Scenario2"],
    narrow_on_page$reservation_wage[narrow_before_page$scenario == "Scenario2"]
  )

  kbl(
    tibble(Scenarios = c("Scenario 1", "Scenario 2"), !!ptestlabel:=c(nvn1, nvn2)),
    "latex",
    escape=FALSE,
    booktabs=T,
    caption=caption,
    align="lr",
    label=label
  ) %>%
    kable_styling(font_size = 12) %>%
    write(filename)
}

tex_for_narrow_v_narrow(
  ptest=wilcox_pvalue,
  ptestlabel="Wilcoxon Test",
  caption="Between-treatment p-values for BOTH when information on baseline is presented for the first time right before the first choice or exactly on the first choice page. Based on two-sided Wilcoxon rank-sum tests, treating each individual in each scenario as a single independent observation.",
  filename="mwu_narrow_v_narrow.tex",
  label="mwu_narrow_v_narrow"
)

tex_for_narrow_v_narrow(
  ptest=t_pvalue,
  ptestlabel="t-test",
  caption="Between-treatment p-values for BOTH when information on baseline is presented for the first time right before the first choice or exactly on the first choice page. Based on two-sided Wilcoxon rank-sum tests, treating each individual in each scenario as a single independent observation.",
  filename="t_test_narrow_v_narrow.tex",
  label="t_test_narrow_v_narrow"
)

## BEFORE and AFTER treatments

means_data(consistent_dfba) %>%
  mutate(treatment = factor(treatment, levels = c("BOTH", "BEFORE", "AFTER", "NONE"), ordered = TRUE)) %>%
  arrange(scenario, treatment) %>%
  select(-scenario) %>%
  table_of_means_data(
    caption="Means of treatments with identical presentation and outcomes but without framing (BOTH), or framing additional tasks as before (BEFORE), or after (AFTER) the baseline tasks, or combining all outcomes (NONE)",
    ntreatments = 4,
    label="before_after"
  ) %>%
  write("before_after.tex")

consistent_dfba_on_page <- consistent_dfba %>%
  filter(!(treatment == "BOTH" & date(session_date) <= january1))

means_data(consistent_dfba_on_page) %>%
  mutate(treatment = factor(treatment, levels = c("BOTH", "BEFORE", "AFTER"), ordered = TRUE)) %>%
  arrange(scenario, treatment) %>%
  select(-scenario) %>%
  table_of_means_data(
    caption= "Restricting BOTH to sessions done with presentation on choice page only. Means of treatments with identical presentation and outcomes but without framing (BOTH), or framing additional tasks as before (BEFORE) or after (AFTER) the baseline tasks",
    ntreatments = 3,
    label="before_after_on_page"
  ) %>%
  write("before_after_on_page.tex")

## Tests with BEFORE AFTER

tex_for_pvalue(
  consistent_dfba,
  wilcox_pvalue,
  caption="Between-treatment p-values for BOTH, BEFORE, AFTER, and NONE treatments based on two-sided Wilcoxon rank-sum tests, treating each individual in each scenario as a single independent observation.",
  label="mwu_before_after",
  filename="mwu_before_after.tex",
  treatments = c("BOTH", "BEFORE", "AFTER", "NONE")
)

tex_for_pvalue(
  consistent_dfba,
  t_pvalue,
  caption="Between-treatment p-values for BOTH, BEFORE, AFTER, and NONE treatments based on two-sided t-tests, treating each individual in each scenario as a single independent observation.",
  label="t_test_before_after",
  filename="t_test_before_after.tex",
  treatments=c("BOTH", "BEFORE", "AFTER", "NONE")
)

tex_for_pvalue(
  consistent_dfba_on_page,
  wilcox_pvalue,
  caption=" Between-treatment p-values for BOTH, BEFORE, AFTER, and NONE treatments based on two-sided Wilcoxon rank-sum tests, treating each individual in each scenario as a single independent observation. Restricted to choices of BOTH when baseline is mentioned on choice page first.",
  label="mwu_before_after_on_page",
  filename="mwu_before_after_on_page.tex",
  treatments = c("BOTH", "BEFORE", "AFTER", "NONE")
)

tex_for_pvalue(
  consistent_dfba_on_page,
  t_pvalue,
  caption="Between-treatment p-values for BOTH, BEFORE, AFTER, and NONE treatments based on two-sided t-tests, treating each individual in each scenario as a single independent observation. Restricted to sessions of BOTH when the endowment is mentioned on choice page first.",
  label="t_test_before_after_on_page",
  filename="t_test_before_after_on_page.tex",
  treatments=c("BOTH", "BEFORE", "AFTER", "NONE")
)

df_balanced <- consistent_df_all %>%
  select(session_date, everything()) %>%
  filter(session_date < dmy("20122019") | ((session_date > dmy("24122019")) & (session_date < dmy("30012020")))) %>%
  rbind(filter(consistent_df_all, treatment %in% c("BEFORE", "AFTER")))

tex_for_pvalue(
  df_balanced,
  wilcox_pvalue,
  caption=" Between-treatment p-values for NONE, BOTH, and MONEY/LOW treatments based on two-sided Wilcoxon rank-sum tests, treating each individual in each scenario as a single independent observation. Restricted to sessions in which these three treatments were balanced.",
  label="mwu_balanced",
  filename="mwu_balanced.tex",
  treatments=main_treatments
)

tex_for_pvalue(
  df_balanced,
  t_pvalue,
  caption=" Between-treatment p-values for NONE, BOTH, and MONEY/LOW treatments based on two-sided t-test, treating each individual in each scenario as a single independent observation. Restricted to sessions in which these three treatments were balanced.",
  label="t_test_balanced",
  filename="t_test_balanced.tex",
  treatments = main_treatments
)

## Money left on the table
## Give the amount of money left on table as difference between BOTH and NONE by Scenario, and by gender, with std errors in parentheses
money_diff <- function(treatments) {
  function(df) {
    t.test(df$reservation_wage [ df$treatment == treatments[1] ], df$reservation_wage [ df$treatment == treatments[2] ])
  }
}

money_diff_low <- function(df) {
  mean(df$reservation_wage [ df$treatment == "NONE" ], na.rm = TRUE) - mean(df$reservation_wage [ df$treatment == "MONEY/LOW" ], na.rm = TRUE)
}

library(broom)

ttests_pooled_gender <- main_df %>%
  nest(-scenario) %>%
  mutate(
    tt = map(data, money_diff(c("NONE", "BOTH"))),
    diff_with_low = map(data, money_diff_low),
    results = map(tt, glance),
    stderr = map(tt, "stderr")
  ) %>%
  unnest(c(results, stderr, diff_with_low)) %>%
  mutate(gender = "Pooled")

ttests_by_gender <- main_df %>%
  filter(gender != "Other") %>%
  nest(-scenario, -gender) %>%
  mutate(
    tt = map(data, money_diff(c("NONE", "BOTH"))),
    diff_with_low = map(data, money_diff_low),
    results = map(tt, glance),
    stderr = map(tt, "stderr")
  ) %>%
  unnest(c(results, stderr, diff_with_low))

# The highest average reservation wage for 15 more tasks is 2.99, which we round to 3.00
# Hence an upper bound for the average cost per task is 3.00/15, which leads to a lower bound for the number of tasks that are equivalent to the money lost.

avg_cost_per_task <- 2.99/15

summary_time_per_task <- consistent_df %>%
  select(time_encryption_main_task, totaltask) %>%
  filter(totaltask > 14) %>%
  mutate(individual_time_per_task = time_encryption_main_task/totaltask) %>%
  summarise(
    avg_time_per_task = mean(individual_time_per_task),
    median_time_per_task = median(individual_time_per_task),
    stderr = sd(individual_time_per_task)
  )

rbind(ttests_pooled_gender, ttests_by_gender) %>%
  select(scenario, gender, estimate, stderr, diff_with_low) %>%
  mutate(
    `Task equivalent`            = estimate/avg_cost_per_task,
    `Time equivalent (in secs)`  = sprintf( (estimate/avg_cost_per_task)*summary_time_per_task$avg_time_per_task, fmt =  '%1.0f')
  ) %>%
  rename(
    Scenario                        = scenario,
    Gender                          = gender,
    `$\\Delta$`                     = estimate,
    `Std.Err.`                      = stderr,
    `$\\hat{\\Delta}$`               = diff_with_low
  ) %>%
  kbl(
    "latex",
    escape=FALSE,
    booktabs=T,
    caption="The table reports $\\Delta$, the change in reservation wages between NONE and BOTH. The highest reservation wage for 15 more tasks is 2.99 across all treatments and scenarios, so that $2.99/15 \\approx 0.20$ is an upper bound for the average cost per task. Using this, we can convert $\\Delta$ into task-equivalents by $\\Delta / 0.20$, and the cost in time-equivalents (in seconds) by $(\\Delta / 0.20) / 46$, since the average time taken for a task is 46 seconds. \\emph{$\\hat{\\Delta}$} stands for $m_{NONE} - m_{MONEY/LOW}$: the change in the marginal disutility of doing 15 extra tasks on top of a low vs on top of a high baseline. Under full narrow bracketing, \\emph{$\\hat{\\Delta}$} and $\\Delta$ should be equal.",
    align="llrrrrr",
    digits=2,
    label="cost_of_bracketing"
  ) %>%
    kable_styling(font_size = 12) %>%
    # FIXME: For final version, don't put gender info in table as well, but for now leave to avoid mistakes in hard-coded names below
    pack_rows("Pooled", 1, 2) %>%
    pack_rows("Female", 3, 4) %>%
    pack_rows("Male", 5, 6) %>%
    write("table_cost_of_bracketing.tex")

## Tables and comparisons for BOTH and MONEY
## CONTINUE NEXT TIME

money_and_none_means <- consistent_df %>%
  filter(treatment %in% c("NONE", "MONEY", "BOTH")) %>%
  nest(-scenario, -treatment) %>%
  mutate(
    Mean = map(data, ~ mean(.$reservation_wage, na.rm = T)),
    N = map(data, ~ dim(.)[[1]])
  ) %>%
  mutate(
    stddev = map(data, ~ sd(.$reservation_wage, na.rm = T))
  ) %>%
  unnest(c(Mean, N, stddev)) %>%
  mutate(
    `Std. Err.` = stddev/sqrt(N),
    Scenario = scenario,
    Treatment = treatment
    ) %>%
  select(Scenario, Treatment, Mean, N, `Std. Err.`) %>%
  kbl(
    "latex",
    escape=FALSE,
    booktabs=T,
    caption="Mean comparisons between NONE, MONEY, and BOTH",
    align="llrrr",
    digits=2,
    label="means_work_bracketing"
  ) %>%
    kable_styling(font_size = 12) %>%
    write("means_work_bracketing.tex")

# Test whether people bracket work broadly

work_diff_low_partial <- function(df) {
  mean(df$reservation_wage [ df$treatment == "MONEY" ], na.rm = TRUE) - mean(df$reservation_wage [ df$treatment == "MONEY/LOW" ], na.rm = TRUE)
}

money_both_table_pooled_gender <- consistent_df %>%
  filter(treatment %in% c("MONEY", "BOTH", "MONEY/LOW")) %>%
  nest(-scenario) %>%
  mutate(
    diff = map(data, money_diff(c("MONEY", "BOTH"))),
    results = map(diff, glance),
    stderr = map(diff, "stderr"),
    diff_with_low_partial = map(data, work_diff_low_partial)
  ) %>%
  unnest(c(results, stderr, diff_with_low_partial)) %>%
  mutate(gender = "Pooled") %>%
  select(scenario, gender, estimate, stderr, diff_with_low_partial, p.value)

money_both_table_by_gender <- consistent_df %>%
  filter(treatment %in% c("MONEY", "BOTH", "MONEY/LOW")) %>%
  # Only 3 observations with gender == "Other", which breaks tests, so drop
  filter(gender != "Other") %>%
  nest(-scenario, -gender) %>%
  mutate(
    diff = map(data, money_diff(c("MONEY", "BOTH"))),
    results = map(diff, glance),
    stderr = map(diff, "stderr"),
    diff_with_low_partial = map(data, work_diff_low_partial)
  ) %>%
  unnest(c(results, stderr, diff_with_low_partial)) %>%
  select(scenario, gender, estimate, stderr, diff_with_low_partial, p.value)

rbind(money_both_table_pooled_gender, money_both_table_by_gender) %>%
  #mutate(
  #  p.value = pvalue_to_char(p.value)
  #) %>%
  select(-p.value) %>%
  rename(
    Scenario = scenario,
    Gender = gender,
    `$\\Delta$` = estimate,
    `Std. Err.` = stderr,
    `$\\hat{\\Delta}$` = diff_with_low_partial
  ) %>%
  kbl(
    "latex",
    escape=FALSE,
    booktabs=T,
    caption="The table reports $\\Delta$, the change in reservation wages between MONEY and BOTH. Under broad bracketing of work, $\\Delta$ should be $0$. \\emph{$\\hat{\\Delta}$} stands for $m_{MONEY} - m_{MONEY/LOW}$. When $\\hat{\\Delta} = 0$, we cannot identify full narrow bracketing from broad bracketing of work.",
    align="llrrr",
    digits=2,
    label="cost_of_work_bracketing"
  ) %>%
    kable_styling(font_size = 12) %>%
    # FIXME: For final version, don't put gender info in table as well, but for now leave to avoid mistakes in hard-coded names below
    pack_rows("Pooled", 1, 2) %>%
    pack_rows("Female", 3, 4) %>%
    pack_rows("Male", 5, 6) %>%
    write("table_cost_of_work_bracketing.tex")

# Test whether people bracket money alone broadly

money_diff_low_partial <- function(df) {
  mean(df$reservation_wage [ df$treatment == "BOTH" ], na.rm = TRUE) - mean(df$reservation_wage [ df$treatment == "MONEY/LOW" ], na.rm = TRUE)
}

money_none_table_pooled_gender <- consistent_df %>%
  filter(treatment %in% c("MONEY", "BOTH", "MONEY/LOW", "NONE")) %>%
  nest(-scenario) %>%
  mutate(
    diff = map(data, money_diff(c("MONEY", "NONE"))),
    results = map(diff, glance),
    stderr = map(diff, "stderr"),
    diff_with_low_partial = map(data, money_diff_low_partial)
  ) %>%
  unnest(c(results, stderr, diff_with_low_partial)) %>%
  mutate(gender = "Pooled") %>%
  select(scenario, gender, estimate, stderr, diff_with_low_partial, p.value)

money_none_table_by_gender <- consistent_df %>%
  filter(treatment %in% c("MONEY", "BOTH", "MONEY/LOW", "NONE")) %>%
  # Only 3 observations with gender == "Other", which breaks tests, so drop
  filter(gender != "Other") %>%
  nest(-scenario, -gender) %>%
  mutate(
    diff = map(data, money_diff(c("MONEY", "NONE"))),
    results = map(diff, glance),
    stderr = map(diff, "stderr"),
    diff_with_low_partial = map(data, money_diff_low_partial)
  ) %>%
  unnest(c(results, stderr, diff_with_low_partial)) %>%
  select(scenario, gender, estimate, stderr, diff_with_low_partial, p.value)

rbind(money_none_table_pooled_gender, money_none_table_by_gender) %>%
  # mutate(
  #   p.value = pvalue_to_char(p.value)
  # ) %>%
  select(-p.value) %>%
  rename(
    Scenario = scenario,
    Gender = gender,
    `$\\Delta$` = estimate,
    `Std. Err.` = stderr,
    `$\\hat{\\Delta}$` = diff_with_low_partial
  ) %>%
  kbl(
    "latex",
    escape=FALSE,
    booktabs=T,
    caption="The table reports $\\Delta$, the change in reservation wages between MONEY and NONE. Broad bracketing of money requires $\\Delta = 0$. \\emph{$\\hat{\\Delta}$} stands for $m_{BOTH} - m_{MONEY/LOW}$. Full narrow bracketing requires $\\hat{\\Delta} = 0$.",
    align="llrrr",
    digits=2,
    label="cost_of_money_bracketing"
  ) %>%
    kable_styling(font_size = 12) %>%
    # FIXME: For final version, don't put gender info in table as well, but for now leave to avoid mistakes in hard-coded names below
    pack_rows("Pooled", 1, 2) %>%
    pack_rows("Female", 3, 4) %>%
    pack_rows("Male", 5, 6) %>%
    write("table_cost_of_money_bracketing.tex")


path_df_NONE_LOWMONEY <- function(sc) {
  offset <- ifelse(sc == "Scenario1", 0, -0.10)
  tibble(
    treatment = c("NONE", "NONE", "MONEY/LOW", "MONEY/LOW"),
    res_wage = c(3.40, 3.50, 3.50, 3.40) - offset,
    scenario = rep(sc, 4)
  )
}

path_df_NONE_BOTH <- function(sc) {
  offset <- ifelse(sc == "Scenario1", 0.10, 0)
  tibble(
    treatment = c("NONE", "NONE", "BOTH", "BOTH"),
    res_wage = c(3.25, 3.35, 3.35, 3.25) - offset,
    scenario = rep(sc, 4)
  )
}

path_df_BOTH_LOWMONEY <- function(sc) {
  offset <- ifelse(sc == "Scenario1", 0.6, 0.25)
  tibble(
    treatment = c("BOTH", "BOTH", "MONEY/LOW", "MONEY/LOW"),
    res_wage = c(3.25, 3.35, 3.35, 3.25) - offset,
    scenario = rep(sc, 4)
  )
}

path_df_BOTH_MONEY <- function(sc) {
  offset <- ifelse(sc == "Scenario1", 0.45, 0.25)
  tibble(
    treatment = c("BOTH", "BOTH", "MONEY", "MONEY"),
    res_wage = c(3.25, 3.35, 3.35, 3.25) - offset,
    scenario = rep(sc, 4)
  )
}

path_df_NONE_MONEY <- function(sc) {
  offset <- ifelse(sc == "Scenario1", 0, -0.10)
  tibble(
    treatment = c("NONE", "NONE", "MONEY", "MONEY"),
    res_wage = c(3.40, 3.50, 3.50, 3.40) - offset,
    scenario = rep(sc, 4)
  )
}


# FIXME: p-values are hardcoded below, that's not ideal
ann_text_main <- tibble(
  treatment = c("NONE", "NONE", "BOTH", "NONE", "NONE", "BOTH"),
  res_wage = c(3.60, 3.35, 2.85, 3.70, 3.45, 3.20),
  lab = c("            < 0.001", "   < 0.001", "    0.106",
          "             0.120", "     0.040", "     0.753"),
  scenario = factor(
    c(rep("Scenario1", 3), rep("Scenario2", 3)),
    c("Scenario1", "Scenario2")
  )
)

ann_text_partial <- tibble(
  treatment = c("NONE", "NONE", "BOTH", "NONE", "NONE", "BOTH"),
  res_wage = c(3.60, 3.35, 3.00, 3.70, 3.45, 3.20),
  lab = c("            0.012", "     < 0.001", "     0.002",
          "            < 0.001", "     0.040", "     0.067"),
  scenario = factor(
    c(rep("Scenario1", 3), rep("Scenario2", 3)),
    c("Scenario1", "Scenario2")
  )
)

scenario1_label <- "Reservation wages in Scenario 1"
scenario2_label <- "Reservation wages in Scenario 2"

scenario1_label_partial <- "Scenario 1: \n d(30) - d(15) in NONE and MONEY"
scenario2_label_partial <- "Scenario 2: \n d(45) - d(30) in NONE and MONEY"

pvalue_plot <- main_df %>%
  means_data() %>%
  ggplot(mapping = aes(x = treatment, y = res_wage)) +
  geom_col(aes(fill = treatment)) +
  geom_errorbar(aes(ymin = res_wage - 2*standard_error, ymax = res_wage + 2*standard_error), width = 0.3) +
  facet_wrap(~scenario, labeller = as_labeller(c(Scenario1 = scenario1_label, Scenario2 = scenario2_label))) +
  geom_path(aes(group = scenario), path_df_NONE_LOWMONEY("Scenario1")) +
  geom_path(aes(group = scenario), path_df_NONE_LOWMONEY("Scenario2")) +
  geom_path(aes(group = scenario), path_df_NONE_BOTH("Scenario1")) +
  geom_path(aes(group = scenario), path_df_NONE_BOTH("Scenario2")) +
  geom_path(aes(group = scenario), path_df_BOTH_LOWMONEY("Scenario1")) +
  geom_path(aes(group = scenario), path_df_BOTH_LOWMONEY("Scenario2")) +
  geom_text(data = ann_text_main, aes(label = lab), hjust=0) +
  labs(x = "Treatment", y = "Reservation Wage ($)", fill = "Treatment")

ggsave("bar_plot_means.png", pvalue_plot)

# FIXME: path_df_BOTH_MONEY and path_df_NONE_MONEY are copy-pasted from money low, so need to be adjusted
# FIXME: ann_text_partial has p-values from MONEY/LOW not MONEY, need to fix

pvalue_plot_partial <- partial_df %>%
  means_data() %>%
  ggplot(mapping = aes(x = treatment, y = res_wage)) +
  geom_col(aes(fill = treatment)) +
  geom_errorbar(aes(ymin = res_wage - 2*standard_error, ymax = res_wage + 2*standard_error), width = 0.3) +
  facet_wrap(~scenario, labeller = as_labeller(c(Scenario1 = scenario1_label_partial, Scenario2 = scenario2_label_partial))) +
  geom_path(aes(group = scenario), path_df_NONE_MONEY("Scenario1")) +
  geom_path(aes(group = scenario), path_df_NONE_MONEY("Scenario2")) +
  geom_path(aes(group = scenario), path_df_NONE_BOTH("Scenario1")) +
  geom_path(aes(group = scenario), path_df_NONE_BOTH("Scenario2")) +
  geom_path(aes(group = scenario), path_df_BOTH_MONEY("Scenario1")) +
  geom_path(aes(group = scenario), path_df_BOTH_MONEY("Scenario2")) +
  geom_text(data = ann_text_partial, aes(label = lab), hjust=0) +
  labs(x = "Treatment", y = "Reservation Wage ($)", fill = "Treatment")

ggsave("bar_plot_means_partial.png", pvalue_plot_partial)
