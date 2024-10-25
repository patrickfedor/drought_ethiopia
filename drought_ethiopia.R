system("gsutil ls gs://cmip6_data/CEDA_CMIP6_Downscaled/subset/Ethiopia/", intern = T)





# SET UP ------------------------------------------------------------

library(tidyverse)
library(stars)
library(furrr)
library(future.apply)

options(future.fork.enable = T)
plan(multicore)
options(future.globals.maxSize = 1000 * 1024^2)

# Set up directories:

dir_gs <- "gs://cmip6_data/CEDA_CMIP6_Downscaled/subset/Ethiopia"

dir_proj <- "/mnt/pers_disk/ethiopia_drought"
fs::dir_create(dir_proj)

dir_data <- str_glue("{dir_proj}/data")
fs::dir_create(dir_data)

dir_res <- str_glue("{dir_proj}/results")
fs::dir_create(dir_res)


# Other global vars:

# time periods for analyses
periods <- 
  list(c(1986,2015),
       c(2040,2080))

# time periods for results
periods_f <- 
  list(c(2041,2060),
       c(2061,2080))




# ASSEMBLE TABLE OF FILES -------------------------------------------

# obtain all file names
ff <- 
  c("pr", "tasmax", "tasmin") %>% 
  map(\(x) {
    system("gsutil ls gs://cmip6_data/CEDA_CMIP6_Downscaled/subset/Ethiopia/", intern = T) %>%
      str_subset(paste0(x, ".*\\.nc$"))
  })


# table of files; details in cols
tb_ff <- 
  tibble(file = unlist(ff) %>% fs::path_file(),
         model = str_split(file, "_", simplify = T) [,4],
         scenario = str_split(file, "_", simplify = T) [,6],
         variable = str_split(file, "_", simplify = T) [,2]) %>% 
  filter(scenario == "ssp585")



# Loop through models
walk(unique(tb_ff$model), \(mod) {

  print(str_glue("PROCESSING MODEL {mod}"))

  
  tb_ff_sub <- tb_ff %>% 
    filter(model == mod)

  # DOWNLOAD DATA ---------------------------------------------------
  
  tb_ff_sub %>%
    future_pwalk(\(file, variable, ...) {
      gs_file_path <- str_glue("{dir_gs}/{file}")
      local_file_path <- str_glue("{dir_data}/{file}")

      print(str_glue("Downloading file from {gs_file_path}"))
    })
      
  # AGGREGATE TO MONTHLY --------------------------------------------
  
  tb_ff_sub$file %>% future_walk(\(f) {
    
    local_file_path <- str_glue("{dir_data}/{f}")
    aggregated_file <- str_replace(f, ".nc", "_mon.nc")

    # Precipitation (pr): total
    # Temperature (tasmax, tasmin): mean
    
    if (str_detect(f, "_pr_")) {
      system(str_glue("cdo monsum {local_file_path} {dir_data}/{aggregated_file}"),
             ignore.stdout = T, ignore.stderr = T)
    } else {
      system(str_glue("cdo monmean {local_file_path} {dir_data}/{aggregated_file}"),
             ignore.stdout = T, ignore.stderr = T)
    }

  })

  # Update table of file names (monthly now)
  tb_ff_sub <- tb_ff_sub %>% 
    mutate(file = str_replace(file, ".nc", "_mon.nc"))

  # SET UP TIME -----------------------------------------------------
  # Extract the time vector
  time_vector <-
    str_glue("{dir_data}/{tb_ff_sub$file[1]}") %>%
    read_ncdf(proxy = T) %>%
    st_get_dimension_values("time") %>%
    as_date()  

  # Index time periods (for ncsub)
  time_ind <- periods %>% map(\(p) {
    c(which(year(time_vector) == p[1]) %>% first(),
      which(year(time_vector) == p[2]) %>% last())
  })

  # IMPORT DATA -----------------------------------------------------

  ss <-
    str_glue("{dir_data}/{tb_ff_sub$file}") %>%
    set_names(tb_ff_sub$variable) %>%
    imap(\(f, i) {

      s <-
        time_ind %>%
        future_map(\(ti) {

          f %>%
            read_ncdf(
              #  only necessary time steps
              ncsub = cbind(start = c(1, 1,  ti[1]),
                            count = c(NA,NA, ti[2]-ti[1]+1)),
              proxy = F) %>%
            suppressMessages() %>%
            setNames("v")

        }) %>%
        do.call(c, .)

      # change tas vars from K to degC
      # precip already in mm
      if (str_detect(i, "tas")) {
        s <-
          s %>%
          mutate(v = v %>% units::set_units(degC))
      }

      return(s)

    })

  
  # CALCULATE PET ---------------------------------------------------
  
  # convert stars to arrays
  a <- ss %>% map(pull)
  
  # Debugging: check the structure of `a` and what variables are present
  print("Debugging: Variables extracted into `a`")
  print(names(a))
  print(dim(a[[1]]))  # Check dimensions for tasmin, tasmax, pr
  
  # extract latitude
  lat <- ss[[1]] %>%
    slice(time, 1) %>%
    stars::st_dim_to_attr(2) %>%
    pull()
  
  # bind latitude at the end of all vars
  # dims will be identical between vars
  # this will enable merging them into a 4-D array
  a <- a %>% map(~abind(.x, lat, along = 3))
  
  # merge into a 4-D array
  a <- do.call(abind, c(a, along = 4)) %>% unname()
  
  # calculate PET using Hargreaves
  s_pet <- a %>%
    future_apply(c(1, 2), \(xx) {
      
      # xx = matrix = timesteps x vars
      print("Debugging: Content of xx before Hargreaves")
      print(xx)
      
      # Check if xx has proper column names (tasmin, tasmax, pr, lat)
      required_vars <- c("tasmin", "tasmax", "pr")
      colnames(xx) <- c("pr", "tasmax", "tasmin", "lat")  # Assign proper column names
      
      # Check if all required variables are present
      if (!all(required_vars %in% colnames(xx))) {
        stop("Missing required variables (tasmin, tasmax, pr) for PET calculation.")
      }
      
      # obtain lat (only length-1 vector needed)
      lat <- xx %>% last() %>% .[1]
      
      # remove lat from the matrix
      xx <- xx[-nrow(xx), ]
      
      # Calculate PET using the Hargreaves method
      PET <- SPEI::hargreaves(
        Tmin = xx[,"tasmin"],
        Tmax = xx[,"tasmax"],
        Pre = xx[,"pr"],
        lat = lat,
        verbose = FALSE,
        na.rm = TRUE
      )
      
      return(PET)
    }) %>%
    aperm(c(2, 3, 1)) 
  
  # convert back to stars object 
  s_pet <- st_as_stars(s_pet) %>%
    setNames("pet")
  
  # homogenize dimensions 
  st_dimensions(s_pet) <- st_dimensions(ss[[1]])
  

})
  
  # # CALCULATE WB ANOMALIES ------------------------------------------
  # 
  # # calculate wb
  # s_wb <- 
  #   c(ss$pr %>% setNames("pr"),
  #     s_pet) %>% 
  #   units::drop_units() %>% 
  #   mutate(wb = pr - pet) %>% 
  #   select(wb)
  # 
  # # index years (to specify baseline)
  # yr_in <- 
  #   s_wb %>% 
  #   st_get_dimension_values("time") %>% 
  #   year() %>% 
  #   unique()
  # 
  # # baseline: 1st period (historical)
  # # 1st year removed: rolling sum
  # 
  # bl_in <- 
  #   (which(yr_in == periods[[1]][1])+1):which(yr_in == periods[[1]][2])
  # 
  # s_wb_anom <- 
  #   s_wb %>% 
  #   st_apply(c(1,2), \(x) {
  #     
  #     if (all(is.na(x))) {
  #       
  #       rep(NA, length(x))
  #       
  #     } else {
  #       
  #       x_cum <- 
  #         x %>% 
  #         zoo::rollsum(k = 12, fill = NA, align = "right")
  #       
  #       m <- 
  #         x_cum %>% 
  #         matrix(ncol = 12, byrow = T)
  #       
  #       m %>%
  #         apply(2, \(xm) round(ecdf(xm[bl_in])(xm), 2)) %>% 
  #         # percentiles based on baseline
  #         t() %>% ##### CHECK !!!!!
  #         as.vector()
  #       
  #     }
  #   },
  #   .fname = "time",
  #   FUTURE = T) %>% 
  #   aperm(c(2,3,1))
  # 
  # # homogenize time dimension
  # st_dimensions(s_wb_anom)[3] <- st_dimensions(s_wb)[3]
  # 
  # # remove 1st year of each period: rolling sum
  # s_wb_anom <- 
  #   s_wb_anom %>% 
  #   filter(!year(time) %in% c(periods[[1]][1], periods[[2]][1]))
  # 
  # 
  # 
  # # CALCULATE PROB OF EXTR DROUGHT ----------------------------------
  # 
  # foo <- 
  #   periods_f %>% # only future time periods
  #   map(\(p) {
  #     
  #     s_wb_anom %>% 
  #       filter(year(time) %in% (p[1]:p[2])) %>% 
  #       st_apply(c(1,2), \(x) {
  #         
  #         # proportion of months under 0.1
  #         # 0.1 = extr. drought perc. threshold
  #         mean(x <= 0.1)
  #         
  #       },
  #       .fname = "prob")
  #     
  #   }) %>% 
  #   do.call(c, .) %>% 
  #   setNames(periods_f %>% 
  #              map(~str_glue("per_{.x[1]}_{.x[2]}")) %>% 
  #              unlist())
  # 
  # # save results
  # #write_rds(foo, str_glue("{dir_res}/{mod}_prob_extr_drought.rds"))
  # 
  # # # delete input data from disk
  # # dir_data %>% 
  # #   fs::dir_ls() %>% 
  # #   fs::file_delete()
  # # 
  # 
})



# CALCULATE FINAL MAPS ----------------------------------------------

s <- 
  dir_res %>% 
  fs::dir_ls() %>% 
  map(read_rds) %>%
  map(merge, name = "period") %>% 
  set_names(unique(tb_ff$model)) %>% 
  {do.call(c, c(., along = "model"))} %>% 
  setNames("prob")


# ensemble means

s_ensmean <- 
  s %>% 
  st_apply(c(1,2,3), \(x) mean(x))

write_stars(s_ensmean, "/mnt/bucket_cmip6/shared_data/Ethiopia/prob_extr_drought_2021_2060.tif")

fs::dir_delete(dir_proj)


# *******************

# coefficient of variation
s_cv <- 
  s %>% 
  st_apply(c(1,2,3), \(x) sd(x)/mean(x)*100) %>% 
  split("period")

# 90% confidence interval
s_ci <- 
  s %>% 
  st_apply(c(1,2,3), \(x) {
    
    if(any(is.na(x))) {
      rep(NA, 2)
    } else {
      
      t.test(x, conf.level = 0.9)$conf.int
      
    }
  },
  .fname = "ci") %>% 
  split("period") %>% 
  aperm(c(2,3,1))


# aggreement 
th <- 0.66 # 2/3
s_aggr <- 
  s %>% 
  st_apply(c(1,2,3), \(x) {
    
    p <- mean(x > 0.1) 
    if_else(p <= (1-th) | p >= th, NA, 1) 
  }
  )


# VISUALIZE

aggr_pol <- 
  bind_rows(
    s_aggr %>%
      split("period") %>% 
      select(1) %>% 
      st_as_sf(as_points = F, merge = T) %>% 
      mutate(period = names(.)[1]) %>% 
      select(-1),
    s_aggr %>%
      split("period") %>% 
      select(2) %>% 
      st_as_sf(as_points = F, merge = T) %>% 
      mutate(period = names(.)[1]) %>% 
      select(-1)
  )


s_ensmean %>% 
  as_tibble() %>% 
  ggplot() +
  geom_raster(aes(lon, lat, fill = prob)) +
  colorspace::scale_fill_continuous_sequential("plasma",
                                               na.value = "transparent",
                                               # breaks = c(0, 0.1, seq(0.15, 0.35, 0.05)),
                                               rev = F,
                                               guide = guide_colorbar(barwidth = 0.6)) +
  geom_sf(data = aggr_pol, fill = "white", alpha = 0.2, color = "grey10") +
  
  coord_sf(expand = F) +
  facet_wrap(~period, ncol = 2) +
  theme(axis.title = element_blank())




