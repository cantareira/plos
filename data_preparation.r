### R code
## Data set preparation and indicators of critical slow-down
### Encoding: UTF-8

###################################################
### R setup
###################################################
require(zoo)
require(plyr)
require(dplyr)
require(caTools)
###################################################
### Reading and transforming data
###################################################
cant <- read.csv("dados.csv", as.is=TRUE) %>%
    filter(manancial=="sistemaCantareira") %>%
        distinct() %>%
            mutate(data = as.Date(data, format="%Y-%m-%d")) %>%
                select(data, volume, pluviometria)
## In Cantareira the usable capacity limit was shift in 18.5% in May 15 2014
## and then in more 10.7% in Oct 23 2014.
cant$volume[cant$data>"2014-05-15"] <- cant$volume[cant$data>"2014-05-15"] - 18.5
cant$volume[cant$data>"2014-10-23"] <- cant$volume[cant$data>"2014-10-23"] - 10.7
### Converting daily time series to zoo object
cant.zoo <- zoo(data.frame(pluv=cant$pluviometria, v.rel=cant$volume,
                           v.abs=(cant$volume+29.2)*9.8155e6,
                           pluv.m=runmean(cant$pluviometria, 30, align="right")), cant$data)
## Inflow and outflow ##
fluxos <- read.csv("data_ocr_cor2.csv", as.is=TRUE) %>%
    mutate(date = as.Date(date, "%Y-%m-%d"),
           inflow=(Jaguari_QNat + Cachoeira_QNat + Atibainha_QNat + PaivaC_QNat)*24*3600,
           outflow=(Jaguari_QJus + Cachoeira_QJus + Atibainha_QJus + PaivaC_QJus +  QESI)*24*3600) %>%
        select(date, inflow, outflow)

fluxos.zoo <- zoo(fluxos[,-1], fluxos$date)
## Merge data on volume and inflow/outflow
cant.plos <- merge(cant.zoo, fluxos.zoo)
## Reading and converting data from the SABESP daily bulletins
## http://site.sabesp.com.br/site/interna/Default.aspx?secaoId=553, available since 2015-01-15
cant.bol <- read.csv("dados_boletins.csv") %>%
    filter(sistema=="Cantareira") %>%
        mutate(data=as.Date(data,"%Y-%m-%d")) %>%
            select(data, vabs)
cant.bol <- zoo(cant.bol$vabs, cant.bol$data)
## Replaces absolute volumes for the more exact figures available since 2015-01-15
cant.plos$v.abs[time(cant.bol)] <- cant.bol*1000000
cant.plos$v.rel[time(cant.bol)] <- cant.plos$v.abs[time(cant.bol)]/9.8155e6 - 29.2
## Final changes
cant.plos <- cant.plos %>% 
    window(start=min(time(fluxos.zoo)), end=max(time(cant.zoo))) %>% ## Croping ends that miss data in one or other series
        na.approx() # Approximating NA values (few)
            
## Selecting data from 2013-01-01
cant.plos.1y <- window(cant.plos, start=max(time(cant.plos))-365)







