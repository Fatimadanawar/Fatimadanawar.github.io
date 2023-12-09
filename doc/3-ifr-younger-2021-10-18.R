## Infection fatality rate of COVID-19 in community-dwelling populations 
## with emphasis on the elderly: An overview 
## Cathrine Axfors, John P.A. Ioannidis
## OSF registration: https://osf.io/47cgb
## medRxiv preprint: https://doi.org/10.1101/2021.07.08.21260210


## README Script 3
# Description: Produces results for main text and Figure 3. Figures were edited for minor 
# layout changes in a separate software.
# Note that "Results in text" section gives line-by-line output, nothing more neat.


# Preparations (hard code) ----
# Set your working directory
# Download and save there "ifr-data3-younger.xlsx"

today <- "2021-10-18"


# Load files ----
agebins <- readxl::read_xlsx("C:/Users/Laptop/Downloads/Workflow-Portfolio/data/ifr-data3-younger.xlsx")

# Load packages ----
library(dplyr)
library(ggplot2)
library(RColorBrewer)


### Data management ----

# infected_agebins	----
#Number of infected people in the age bin for date 1.
#Estimated by multiplying the adjusted estimate (and if unavailable, the unadjusted) of
#seroprevalence in the age bin with pop_agebins.

agebins <- purrr::modify_at(agebins, c("adj_sero_agebins", "crude_sero_agebins",
                                       "sero_agebins_g1", "sero_agebins_g2", "sero_agebins_g3",
                                       "pop_agebins_g1", "pop_agebins_g2", "pop_agebins_g3",
                                       "sero_agebins_g4", "pop_agebins_g4",
                                       "pop_agebins", "number_sampled"), 
                            ~as.numeric(.x))

agebins <- mutate(agebins, infected_agebins = case_when(
  !is.na(adj_sero_agebins) ~ adj_sero_agebins/100*pop_agebins,
  !is.na(sero_agebins_g4) ~
    sero_agebins_g4/100*pop_agebins_g4 +
    sero_agebins_g3/100*pop_agebins_g3 + 
    sero_agebins_g2/100*pop_agebins_g2 + 
    sero_agebins_g1/100*pop_agebins_g1,
  !is.na(sero_agebins_g3) ~ 
    sero_agebins_g3/100*pop_agebins_g3 + 
    sero_agebins_g2/100*pop_agebins_g2 + 
    sero_agebins_g1/100*pop_agebins_g1,
  !is.na(sero_agebins_g2) ~
    sero_agebins_g2/100*pop_agebins_g2 + 
    sero_agebins_g1/100*pop_agebins_g1,
  !is.na(crude_sero_agebins) ~ crude_sero_agebins/100*pop_agebins,
  TRUE ~ NA_real_
))


# deaths_agebins_date1 ----
#COVID-19 deaths in age bin for date 1, estimated using deaths_agebins_date2, age_date2,
#and deaths_date1

agebins <- mutate(agebins, deaths_agebins_date1 = case_when(
  age_date2 == "Not applicable" ~ as.numeric(deaths_agebins_date2),
  !is.na(as.numeric(deaths_agebins_date2)) ~ 
    as.numeric(deaths_agebins_date2)/as.numeric(age_date2)*as.numeric(deaths_date1),
  TRUE ~ NA_real_
))


# ifr_uncorr_agebins ----
#IFR in age bins, uncorrected

agebins <- mutate(agebins, ifr_uncorr_agebins = case_when(
  infected_agebins == 0 ~ 0,
  !is.na(deaths_agebins_date1) ~ deaths_agebins_date1/infected_agebins,
  TRUE ~ NA_real_
))


# ifr_corr_agebins ----
#IFR in age bins, corrected

agebins <- mutate(agebins, antibody_type2 = case_when(
  grepl("total|pan.Ig|IgG, IgM, IgA", antibody_type, ignore.case = T) ~ "IgG/IgM/IgA",
  grepl("IgG and.or IgM|IgG, IgM", antibody_type, ignore.case = T) ~ "IgG/IgM",
  grepl("missing|unclear", antibody_type, ignore.case = T) ~ "Missing/Unclear",
  TRUE ~ antibody_type
))

agebins <- mutate(agebins, ifr_corr_agebins = case_when(
  grepl("IgG.IgM.IgA|missing", antibody_type2, ignore.case = T) ~ 
    ifr_uncorr_agebins,
  grepl("^IgG$", antibody_type2, ignore.case = T) ~ ifr_uncorr_agebins/1.1/1.1,
  grepl("IgG.IgM", antibody_type2, ignore.case = T) ~ ifr_uncorr_agebins/1.1
))



# Mid-point of age bin ----
#For Figure 3

agebins <- purrr::modify_at(agebins, c("upper_age", "lower_age"), 
                            ~as.numeric(.x))

agebins$midpoint <- round(agebins$lower_age + 
                            (agebins$upper_age-agebins$lower_age)/2, digits = 1)

#Delete age bins that are not eligible
agebins$width <- agebins$upper_age-agebins$lower_age
agebins <- mutate(agebins, included_plot = case_when(
  grepl("Not include", comment, ignore.case = T) ~ "No",
  width < 20 ~ "Yes",
  TRUE ~ "No"
))

agebins <- subset(agebins, included_plot == "Yes")


# Calculation of medians ----

# location2
agebins <- mutate(agebins, location2 = case_when(
  grepl("England|UK|Scotland", location) ~ "UK",
  grepl("Canada", location) ~ "Canada",
  grepl("India", location) ~ "India",
  grepl("Qatar", location) ~ "Qatar",
  grepl("France", location) ~ "France",
  TRUE ~ location
))


#1. Separate datasets per age bin

agebins0_19 <- filter(agebins, midpoint<20)
agebins20_29 <- filter(agebins, midpoint>19.9 & midpoint<30)
agebins30_39 <- filter(agebins, midpoint>29.9 & midpoint<40)
agebins40_49 <- filter(agebins, midpoint>39.9 & midpoint<50)
agebins50_59 <- filter(agebins, midpoint>49.9 & midpoint<60)
agebins60_69 <- filter(agebins, midpoint>59.9 & midpoint<70)
#test <- agebins0_19[c("lower_age", "upper_age", "midpoint")]
#test <- agebins20_29[c("lower_age", "upper_age", "midpoint")]
#test <- agebins30_39[c("lower_age", "upper_age", "midpoint")]
#test <- agebins40_49[c("lower_age", "upper_age", "midpoint")]
#test <- agebins50_59[c("lower_age", "upper_age", "midpoint")]
#test <- agebins60_69[c("lower_age", "upper_age", "midpoint")]


#2. Sample size weighted averages for same country and age bin

agebins0_19 <- group_by(agebins0_19, location2)
n_groups(agebins0_19)
agebins0_19 <- mutate(agebins0_19, weight = case_when(
  !is.na(number_sampled) ~ number_sampled/sum(number_sampled),
  TRUE ~ 1
))
agebins0_19 <- mutate(agebins0_19, step = ifr_corr_agebins*weight)
agebins0_19 <- mutate(agebins0_19, ifr_country_weighted = sum(step))
agebins0_19 <- ungroup(agebins0_19)
agebins0_19$newmidpoint <- "10"
#test <- agebins0_19[c("location2", "ifr_corr_agebins", "ifr_country_median", "number_sampled", "weight", "ifr_country_weighted")]


agebins20_29 <- group_by(agebins20_29, location2)
n_groups(agebins20_29)
agebins20_29 <- mutate(agebins20_29, weight = case_when(
  !is.na(number_sampled) ~ number_sampled/sum(number_sampled),
  TRUE ~ 1
))
agebins20_29 <- mutate(agebins20_29, step = ifr_corr_agebins*weight)
agebins20_29 <- mutate(agebins20_29, ifr_country_weighted = sum(step))
agebins20_29 <- ungroup(agebins20_29)
agebins20_29$newmidpoint <- "25"
#test <- agebins20_29[c("location2", "ifr_corr_agebins", "ifr_country_median", "number_sampled", "weight", "ifr_country_weighted")]


agebins30_39 <- group_by(agebins30_39, location2)
n_groups(agebins30_39)
agebins30_39 <- mutate(agebins30_39, weight = case_when(
  !is.na(number_sampled) ~ number_sampled/sum(number_sampled),
  TRUE ~ 1
))
agebins30_39 <- mutate(agebins30_39, step = ifr_corr_agebins*weight)
agebins30_39 <- mutate(agebins30_39, ifr_country_weighted = sum(step))
agebins30_39 <- ungroup(agebins30_39)
agebins30_39$newmidpoint <- "35"
#test <- agebins30_39[c("location2", "ifr_corr_agebins", "ifr_country_median", "number_sampled", "weight", "ifr_country_weighted")]


agebins40_49 <- group_by(agebins40_49, location2)
n_groups(agebins40_49)
agebins40_49 <- mutate(agebins40_49, weight = case_when(
  !is.na(number_sampled) ~ number_sampled/sum(number_sampled),
  TRUE ~ 1
))
agebins40_49 <- mutate(agebins40_49, step = ifr_corr_agebins*weight)
agebins40_49 <- mutate(agebins40_49, ifr_country_weighted = sum(step))
agebins40_49 <- ungroup(agebins40_49)
agebins40_49$newmidpoint <- "45"
#test <- agebins40_49[c("location2", "ifr_corr_agebins", "ifr_country_median", "number_sampled", "weight", "ifr_country_weighted")]


agebins50_59 <- group_by(agebins50_59, location2)
n_groups(agebins50_59)
agebins50_59 <- mutate(agebins50_59, weight = case_when(
  !is.na(number_sampled) ~ number_sampled/sum(number_sampled),
  TRUE ~ 1
))
agebins50_59 <- mutate(agebins50_59, step = ifr_corr_agebins*weight)
agebins50_59 <- mutate(agebins50_59, ifr_country_weighted = sum(step))
agebins50_59 <- ungroup(agebins50_59)
agebins50_59$newmidpoint <- "55"
#test <- agebins50_59[c("location2", "ifr_corr_agebins", "ifr_country_median", "number_sampled", "weight", "ifr_country_weighted")]


agebins60_69 <- group_by(agebins60_69, location2)
n_groups(agebins60_69)
agebins60_69 <- mutate(agebins60_69, weight = case_when(
  !is.na(number_sampled) ~ number_sampled/sum(number_sampled),
  TRUE ~ 1
))
agebins60_69 <- mutate(agebins60_69, step = ifr_corr_agebins*weight)
agebins60_69 <- mutate(agebins60_69, ifr_country_weighted = sum(step))
agebins60_69 <- ungroup(agebins60_69)
agebins60_69$newmidpoint <- "65"
#test <- agebins60_69[c("location2", "ifr_corr_agebins", "ifr_country_median", "number_sampled", "weight", "ifr_country_weighted")]


#3. Add to main dataset

agebins_new <- rbind(agebins0_19, agebins20_29, agebins30_39, agebins40_49, agebins50_59, agebins60_69)
agebins_new <- group_by(agebins_new, newmidpoint)
agebins_new <- filter(agebins_new, !duplicated(location2))
agebins_new <- mutate(agebins_new, ifr_mdn_weighted = median(ifr_country_weighted))
n_groups(agebins_new)
agebins_new <- ungroup(agebins_new)
test <- select(agebins_new, c("newmidpoint", "ifr_mdn_weighted"))
test$ifr_mdn_weighted <- test$ifr_mdn_weighted*100
test <- subset(test, !duplicated(ifr_mdn_weighted))


### Output ----

# "IFR in younger age-strata" ----

length(unique(agebins$study))
length(unique(agebins$location2))

test #Median IFR in the different age groups (they are here defined by their midpoints)

length(unique(agebins0_19$location2))
length(unique(agebins20_29$location2))
length(unique(agebins30_39$location2))
length(unique(agebins40_49$location2))
length(unique(agebins50_59$location2))
length(unique(agebins60_69$location2))


# Output for the purposes of proofreading ----

writexl::write_xlsx(agebins_new, paste0("ifr-younger-all-variables-", today, ".xlsx"))



# Figure 3 ----

agebins_new <- group_by(agebins_new, newmidpoint)
n_groups(agebins_new)
agebins_new <- filter(agebins_new, !duplicated(location2))
agebins_new <- ungroup(agebins_new)

ageplot1 <- ggplot(agebins_new, aes(x = newmidpoint, y = ifr_country_weighted*100))+
  geom_point(aes(color=location2), size=3)+
  theme_bw()+
  geom_line(aes(group=location2, color=location2))+
  theme(panel.grid = element_blank())+
  labs(y = "IFR (%)", x = "Age", color="Country")+
  scale_color_manual(values = c(brewer.pal(9, "Set1")[-6], 
                                rev(brewer.pal(8, "Dark2")),
                                rev(brewer.pal(8, "Accent"))))+
  scale_y_log10()+
  theme(strip.background = element_blank())
pdf(paste0("Figure3-", today, ".pdf"), height = 5.5, width = 7)
plot(ageplot1)
dev.off()

