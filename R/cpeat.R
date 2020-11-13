clymo <-function(data, a = 0.01, b = 0.01, grouping_var = site){
  require(forestmangr)
  require(tidyverse)
  options(warn = -1)
  a <- a
  b <- b
  grouping_var <- enquo(grouping_var)

  data <- data %>%
    group_by(!!grouping_var) %>%
    mutate(cumulative_mass = cumsum((depth - lag(depth, default = 0)) * bulk_density)) %>%
    ungroup()

  nls_table(data,
            cumulative_mass ~ (a / b) * (1 - exp(-b * age)),
            mod_start = c(a = a, b = b),
            output = "table",
            .groups = "site"
  ) %>%
    as.data.frame%>%
    rename(addition =b0, decomp = b1)}


prep <- function(data, data1, grouping_var = site){
  options(warn = -1)
  require(tidyverse)
  require(imputeTS)

  grouping_var <-  enquo(grouping_var)

  data <- data %>%
    group_by(!!grouping_var) %>%
    mutate(t = row_number())%>%
    ungroup()


  data1 <- data1 %>%
    group_by(!!grouping_var) %>%
    mutate(t = row_number())%>%
    ungroup()


  df <- left_join(data, data1) %>%
    na_replace(0)}


acrotelm <- function(data, grouping_var = site) {
  require(tidyverse)
  require(imputeTS)
  options(warn = -1)
  grouping_var <- enquo(grouping_var)
  acrotelm <- data %>%
    group_by(!!grouping_var) %>%
    arrange(!!grouping_var, depth) %>%
    mutate(t = row_number(),
           annual_layers = (max(addition) / max(decomp)) * log((1 + max(decomp) * lead(t)) /
                                                                 (1 + max(decomp) * t)),
           mass_loss = max(decomp) * ((1 / (1 + max(decomp) * lead(t)) - (1 / (1 +max(decomp) * t
           ))) / log((1 + max(decomp) * lead(t)) / (1 + max(decomp) * t))),
           decomp_rate = max(decomp) / (1 + max(decomp) * t),
           peat_decomp = max(decomp) / (1 + max(decomp) * t),
           annual_layers1 = (max(addition) / min(decomp)) * log((1 + min(decomp) * lead(t)) /
                                                                  (1 + min(decomp) * t)),
           mass_loss1 = min(decomp) * ((1 / (1 + min(decomp) * lead(t)) - (1 / (1 +min(decomp) * t
           ))) / log((1 + min(decomp) * lead(t)) / (1 + max(decomp) * t))),
           decomp_rate1 = min(decomp) / (1 + min(decomp) * t),
           peat_decomp1 = min(decomp) / (1 + min(decomp) * t))%>%
    na_replace(0) %>%
    mutate(annual_layers = ifelse(annual_layers == 0, annual_layers1, annual_layers),
           mass_loss = ifelse(mass_loss == 0, mass_loss1, mass_loss),
           decomp_rate = ifelse(decomp_rate == 0, decomp_rate1, decomp_rate),
           peat_decomp = ifelse(peat_decomp == 0, peat_decomp1, peat_decomp)) %>%
    select(!c(annual_layers1, mass_loss1, decomp_rate1, peat_decomp1, depth, age,
              zeroed_age, bulk_density, t, carbon, addition, decomp)) %>%
    ungroup()}


mega_bog <- function(data, bin_size = 100, grouping_var = site) {
  options(warn = -1)
  require(tidyverse)
  require(imputeTS)
  bin_size <- bin_size
  grouping_var <- enquo(grouping_var)
  data1 <- data %>%
    group_by(!!grouping_var) %>%
    arrange(!!grouping_var, depth) %>%
    mutate(time = lead(zeroed_age) - zeroed_age) %>%
    na_replace(0) %>%
    mutate(cumulative_time = cumsum(time),
           cumulative_mass = cumsum((depth - lag(depth, default = 0))
                                    * bulk_density)) %>%
    na_replace(0) %>%
    mutate(
      time_bin = ceiling(cumulative_time / bin_size) * bin_size,
      time_bin = ifelse(time_bin == 0, max(time_bin) + bin_size, time_bin)
    ) %>%
    complete(time_bin = seq(min(time_bin), max(time_bin), by = bin_size)) %>%
    na_interpolation() %>%
    mutate(
      time = lead(zeroed_age) - zeroed_age,
      cumulative_time = cumsum(time),
      RERCA = ((lead(depth) - depth) / (lead(zeroed_age) - zeroed_age)) * bulk_density * carbon * 10000,
      rerca_net_carbon_pool = RERCA * (lead(zeroed_age) - zeroed_age)
    ) %>%
    na_replace(0)

  mega_long <<- data1 %>%
    mutate(
      cum_c_pool_Bottom_Up = rev(cumsum(rev(rerca_net_carbon_pool))),
      cum_c_pool_Top_Down = cum_c_pool_Bottom_Up[1] - cum_c_pool_Bottom_Up,
      dif = time_bin - cumulative_time,
      ndif = time_bin - lead(cumulative_time),
      ndif = ndif * -1,
      ndif1 = time - lag(dif),
      ndif1 = ifelse(ndif1 > 0, ndif1, 0)
    ) %>%
    na_replace(0) %>%
    group_by(time_bin, !!grouping_var) %>%
    mutate(ndif = case_when(ndif != min(ndif) ~ 0, TRUE ~ ndif),
           dif = case_when(dif != min(dif) ~ 0, TRUE ~ dif)) %>%
    ungroup(time_bin) %>%
    mutate(ndif = lag(ndif)) %>%
    na_replace(0) %>%
    mutate(
      time1 = ifelse(ndif > 0, ndif, time),
      time1 = ifelse(dif > 0, dif, time1),
      RERCA1 = lead(RERCA),
      dif = case_when(dif != min(dif) ~ 0, TRUE ~ dif)
    ) %>%
    group_by(!!grouping_var, time_bin) %>%
    do(add_row(., .before = max())) %>%
    ungroup() %>%
    fill(!!grouping_var) %>%
    group_by(!!grouping_var) %>%
    na_replace(0) %>%
    mutate(
      time_bin = ifelse(time_bin == 0, lag(time_bin), time_bin),
      RERCA2 = lead(RERCA),
      time2 = lag(time1),
      RERCA = ifelse(RERCA == 0, RERCA2, RERCA),
      time = ifelse(time == 0, time2, time),
      time = ifelse(ndif1 == 0, time, ndif1),
      time = ifelse(time > 100, time - 100, time)
    ) %>%
    select(!c(dif, ndif, ndif1, time1, RERCA1, RERCA2, time2)) %>%
    ungroup() %>%
    slice(-1)


  data2 <- mega_long %>%
    group_by(time_bin, !!grouping_var) %>%
    mutate(rerca = sum(RERCA * time)) %>%
    slice(1) %>%
    ungroup(time_bin) %>%
    mutate(RERCA_Mean = rerca / bin_size) %>%
    na_replace(0) %>%
    mutate(c_pool = rev(cumsum(rev(RERCA_Mean)))) %>%
    na_replace(0) %>%
    slice(1:(n() - 1)) %>%
    mutate(
      c_pool_clymo = c_pool[1] - c_pool,
      net_carbon_pool = lead(c_pool_clymo) - c_pool_clymo) %>%
    na_replace(0) %>%
    mutate(net_c_release = pmap_dbl(list(seq_len(n()), exp(-0.0008 * time_bin),
                                         exp(-0.0008 * c(0, time_bin[-n()]))),
                                    function(i, tb, tblag)
                                      sum(net_carbon_pool[i:n()] / tb - net_carbon_pool[i:n()] / tblag))
    )

  data3 <- data2 %>%
    group_by(!!grouping_var) %>%
    mutate(time_bin1 = lag(time_bin),
           net_c_uptake = net_carbon_pool / (exp(-max(decomp) * time_bin))) %>%
    na_replace(0) %>%
    mutate(
      net_c_uptake = ifelse(net_c_uptake == 0, net_carbon_pool / (exp(-min(decomp) * time_bin)), net_c_uptake),
      c_relase = net_c_uptake - net_carbon_pool,
      c_release_100 = net_carbon_pool / (exp(-max(decomp) * time_bin)) - net_carbon_pool /
        (exp(-max(decomp) * time_bin1))
    ) %>%
    na_replace(0) %>%
    mutate(
      c_release_100 = ifelse(c_release_100 == 0, net_carbon_pool / (exp(
        -min(decomp) * time_bin
      )) - net_carbon_pool /
        (exp(
          -min(decomp) * time_bin1
        )), c_release_100),
      cum_c_uptake = rev(cumsum(rev(net_c_uptake)))
    ) %>%
    na_replace(0) %>%
    mutate(
      neg_net_c_release = net_c_release * -1,
      cum_c_release = rev(cumsum(rev(net_c_release))),
      neg_cum_c_release = cum_c_release * -1,
      net_c_balance = net_c_uptake - net_c_release,
      cum_c_pool = cum_c_uptake - cum_c_release,
      cum_cohort_c_release = net_c_uptake - net_carbon_pool,
      neg_cum_cohort_c_release = cum_cohort_c_release * -1
    ) %>%
    select(!time_bin1)
}



