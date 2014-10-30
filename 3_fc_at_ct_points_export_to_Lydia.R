library(dplyr)
library(reshape2)
library(ggplot2)

load("ct_fc.RData")
load("ct_periods.RData")

# Convert from n_pixels in traps_fc into hectares
pix_cols <- grepl('^((fc_2001)|(gain)|(lossgain)|(l20))', names(traps_fc))
traps_fc[pix_cols] <- traps_fc[pix_cols]*30^2/(100*100)

ct_years <- group_by(ct_periods, sitecode) %>%
    summarize(ct_start_year=min(year)) %>%
    mutate(fc_start_year=ct_start_year-5,
           fc_end_year=ct_start_year)

traps_loss_lng <- melt(select(traps_fc, sitecode, trap_ID, buffer_m, 
                              starts_with('l20')), variable.name="year", 
                       value.name="loss")
traps_loss_lng$year <- as.numeric(gsub('l', '', traps_loss_lng$year))

# Calculate forest cover loss in 5 year period prior to start of ct deployment
calc_pre_ct_loss <- function(sitecode, year, loss) {
    if (any(is.na(loss))) return(NA)
    start_year <- ct_years$fc_start_year[ct_years$sitecode == sitecode]
    end_year <- ct_years$fc_end_year[ct_years$sitecode == sitecode]
    sum(loss[(year > start_year) & (year <= end_year)])
}
# Calculate forest cover 5 years prior to start of ct deployment
calc_pre_ct_cover <- function(sitecode, trap_ID, buffer_m, year, loss) {
    if (any(is.na(loss))) return(NA)
    fc_2001 <- traps_fc$fc_2001[traps_fc$sitecode == sitecode & 
                                traps_fc$trap_ID == trap_ID &
                                traps_fc$buffer_m == buffer_m]
    start_year <- ct_years$fc_start_year[ct_years$sitecode == sitecode]
    # sum the loss up until the start of 5 year period prior to ct, and 
    # subtract this from 2001 forest cover
    fc_2001 - sum(loss[year <= start_year])
}
ct_loss <- group_by(traps_loss_lng, sitecode, trap_ID, buffer_m) %>%
    summarize(fc_5yr_prior_ct=calc_pre_ct_cover(sitecode[1], trap_ID[1], buffer_m[1], year, loss),
              loss_5yr_period=calc_pre_ct_loss(sitecode[1], year, loss),
              loss_2000_2012=sum(loss, na.rm=TRUE))

save(ct_loss, file="ct_fc_loss.RData")

table(ct_loss$loss_5yr_period)
table(ct_loss$loss_5yr_period > 0, ct_loss$sitecode)

ggplot(ct_loss) +
    geom_bar(aes(loss_5yr_period)) +
    facet_grid(buffer_m~., scales="free")

ggplot(ct_loss) +
    geom_bar(aes(loss_2000_2012)) +
    facet_grid(buffer_m~., scales="free")

ggplot(ct_loss) +
    geom_bar(aes(fc_5yr_prior_ct)) +
    facet_grid(buffer_m~., scales="free")

ggplot(traps_fc) +
    geom_bar(aes(fc_2001)) +
    facet_grid(buffer_m~., scales="free")

table(ct_loss$loss_5yr_period > 0, ct_loss$sitecode)
sum(ct_loss$loss_5yr_period > 0, na.rm=TRUE)
table(ct_loss$loss_2000_2012 > 0, ct_loss$sitecode)
sum(ct_loss$loss_2000_2012 > 0, na.rm=TRUE)

ggplot(traps_loss_lng) +
    geom_bar(aes(factor(year), loss, fill=sitecode), position="dodge", stat="identity") +
    facet_grid(buffer_m~., scales="free")
