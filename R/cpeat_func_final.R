

clymo <-function(data, a, b){
  require(forestmangr)
  require(tidyverse)
  a <- a
  b <- b

  data <- data %>%
    group_by(site) %>%
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


prep <- function(data, data1){
  require(tidyverse)
  require(imputeTS)

  data <- data %>%
  group_by(site) %>%
  mutate(t = row_number())%>%
  ungroup()


data1 <- data1 %>%
  group_by(site) %>%
  mutate(t = row_number())%>%
  ungroup()


df <- left_join(data, data1) %>%
  na_replace(0)}


acrotelm <- function(data) {
  require(tidyverse)
  require(imputeTS)
  acrotelm <- data %>%
    group_by(site) %>%
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


mega_bog <- function(data) {
  options(scipen = 10)
  require(tidyverse)
  require(imputeTS)
  data1 <- data %>%
    group_by(site)%>%
    arrange(site, depth) %>%
    mutate(time = lead(age) - age) %>%
    replace(is.na(.), 0)%>%
    group_by(site) %>%
    mutate(cumulative.time = cumsum(time)) %>%
    na_replace(0) %>%
    mutate(time_bin = ceiling(cumulative.time / 100) * 100,
           time_bin = ifelse(time_bin == 0, max(time_bin) + 100, time_bin))%>%
    complete(time_bin = seq(min(time_bin), max(time_bin), by = 100))%>%
    na_interpolation() %>%
    mutate(time = lead(age) - age) %>%
    replace(is.na(.), 0)%>%
    group_by(site) %>%
    mutate(cumulative.time = cumsum(time),
      PCAR = ((
        lead(depth) - depth
      ) / (lead(age) - age)) * bulk_density * carbon * 10000,
      PCA_NCP = PCAR * (lead(age) - age)
    ) %>%
    na_replace(0)

mega_long <<- data1 %>%
  mutate(CCP_Bottom_Up = rev(cumsum(rev(PCA_NCP))),
                     CCP_Top_Down = CCP_Bottom_Up[1] - CCP_Bottom_Up) %>%
    na_replace(0) %>%
    arrange(site, depth) %>%
     na_replace(0) %>%
    mutate(
      dif = time_bin - cumulative.time,
      ndif = time_bin - lead(cumulative.time),
      ndif = ndif * -1,
      ndif1 = time - lag(dif),
      ndif1 = ifelse(ndif1 > 0, ndif1, 0)
    ) %>%
    group_by(time_bin, site) %>%
    mutate(ndif = case_when(ndif != min(ndif) ~ 0, TRUE ~ ndif),
           dif = case_when(dif != min(dif) ~ 0, TRUE ~ dif))%>%
    ungroup() %>%
    group_by(site) %>%
    mutate(ndif = lag(ndif)) %>%
    na_replace(0) %>%
    mutate(time1 = ifelse(ndif > 0, ndif, time),
      time1 = ifelse(dif > 0, dif, time1),
      PCAR1 = lead(PCAR),
      dif = case_when(dif != min(dif) ~ 0, TRUE ~ dif)) %>%
    group_by(site, time_bin)%>%
    do(add_row(., .before = max())) %>%
  ungroup() %>%
    fill(site) %>%
    group_by(site) %>%
    na_replace(0) %>%
    mutate(time_bin = ifelse(time_bin == 0, lag(time_bin), time_bin))%>%
    mutate(PCAR2 = lead(PCAR),
           time2 = lag(time1),
           PCAR = ifelse(PCAR == 0, PCAR2, PCAR),
           time = ifelse(time == 0, time2, time),
           time = ifelse(ndif1 == 0, time, ndif1),
           time = ifelse(time > 100, time - 100, time)) %>%
  select(!c(dif, ndif, ndif1, time1, PCAR1,PCAR2, time2)) %>%
  ungroup() %>%
  slice(-1)


data2 <- mega_long%>%
  group_by(time_bin, site)%>%
  mutate(PCA = sum(PCAR*time)) %>%
  slice(1) %>%
    ungroup %>%
    group_by(site) %>%
    mutate(
    PCAR_Mean = PCA / 100) %>%
    na_replace(0)%>%
    mutate(cpool = rev(cumsum(rev(PCAR_Mean)))) %>%
    na_replace(0) %>%
    slice(1:(n()-1)) %>%
   mutate(cpool_clymo = cpool[1] - cpool,
    NCP = lead(cpool_clymo) - cpool_clymo) %>%
  na_replace()%>%
    ungroup() %>%
    group_by(site) %>%
    mutate(NCR = pmap_dbl(list(
      seq_len(n()), exp(-0.0008 * time_bin),
      exp(-0.0008 * c(0, time_bin[-n()]))
    ),
    function(i, tb, tblag)
      sum(NCP[i:n()] / tb - NCP[i:n()] / tblag)))

data3 <- data2 %>%
  group_by(site) %>%
  mutate(
    time_bin1 = lag(time_bin)) %>%
na_replace(0) %>%
  mutate(NCU = NCP / (exp(-max(decomp) * time_bin)))%>%
  na_replace(0) %>%
  mutate(NCU= ifelse(NCU == 0, NCP / (exp(-min(decomp) * time_bin)), NCU),
  C_relase = NCU - NCP,
C_release_100 = NCP/(exp(-max(decomp)*time_bin))-NCP/
  (exp(-max(decomp) * time_bin1))) %>%
  na_replace(0) %>%
  mutate(C_release_100 = ifelse(C_release_100 == 0, NCP/(exp(-min(decomp)*time_bin))-NCP/
                                  (exp(-min(decomp) * time_bin1)), C_release_100),
 CCU = rev(cumsum(rev(NCU)))) %>%
  na_replace(0) %>%
  mutate(
    negNCR = NCR * -1,
    CCR = rev(cumsum(rev(NCR))),
    negCCR = CCR * -1,
    NCB = NCU - NCR,
    CCP = CCU - CCR,
    CCCR = NCU - NCP,
    negCCCR = CCCR * -1)}

