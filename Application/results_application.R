############################################################
# Required packages
############################################################

library(loo)
library(data.table)
library(dplyr)
library(igraph)

# Working directory containing the fitted models and application files
setwd(
  "C:/Users/alejo/Dropbox/postdoc20242/random_graphs_new/scripts_random_graphs_prof_mauricio_ordenados/new_simulations/Codigo github biom/Application"
)


############################################################
# Load fitted models and prediction data
############################################################

# Model fitted using a randomly selected spatial graph
res_rand = readRDS(
  "modelo_grafo_2_intercepto_selec_random_2.rds"
)

# Model without spatial dependence, retaining the temporal component
res_time = readRDS(
  "modelo_grafo_2_intercepto_selec_ind_temp_2.rds"
)

# Model fitted using the reported/original spatial graph
res_rep = readRDS(
  "modelo_grafo_2_intercepto_selec_reported_2.rds"
)


# Observations corresponding to the forecasting period
datapred = readRDS(
  "datapred.rds"
)


############################################################
# PSIS-LOO evaluation of candidate spatial graphs
############################################################

### Joining the results of the selected models
listaux = list()
for(i in 1:10){
  outcome3 =  readRDS(paste0("modelo_grafo_2_selec_2_",i,".rds"))
  listaux[[i]] = list(pred_dist = outcome3$pred_dist, log_lik = outcome3$log_lik)
}

h = listaux

mse = 0
elpdind = 0

# Stores the loo object for each candidate graph.
# These objects are later used to estimate model-averaging weights.
arg_weights = list()

se = 0


for (k in 1:length(h)) {
  
  # Compute PSIS-LOO for candidate graph k
  arg_weights[[k]] =
    loo(h[[k]]$log_lik)
  
  # Expected log predictive density (ELPD)
  elpdind[k] =
    arg_weights[[k]]$estimates[
      "elpd_loo",
      "Estimate"
    ]
  
  # Standard error associated with the ELPD estimate
  se[k] =
    arg_weights[[k]]$estimates[
      "elpd_loo",
      "SE"
    ]
}


############################################################
# Rank candidate graphs according to ELPD
############################################################

# Larger ELPD values indicate better predictive performance.
ord <- order(
  elpdind,
  decreasing = TRUE
)

elpdindord <- elpdind[ord]
seord <- se[ord]


############################################################
# Compute model-averaging weights
############################################################

# Stacking weights are estimated from the PSIS-LOO objects.
model_weights_1 =
  loo_model_weights(
    arg_weights
  )

# Pseudo-BMA weights provide an alternative model-averaging scheme.
model_weights_2 =
  loo_model_weights(
    arg_weights,
    method = "pseudobma"
  )


weights =
  as.numeric(model_weights_1)

weights_2 =
  as.numeric(model_weights_2)


############################################################
# Construct model-averaged predictive distributions
############################################################

# Stacking predictive distribution
ypreddistfin = 0

# Pseudo-BMA predictive distribution
ypreddistfin3 = 0

# Weighted log-likelihood matrices
elpd1aux = 0
elpd3aux = 0


for (s in 1:length(weights)) {
  
  # Weighted posterior predictive draws using stacking
  ypreddistfin =
    ypreddistfin +
    weights[s] *
    h[[s]]$pred_dist
  
  # Weighted posterior predictive draws using pseudo-BMA
  ypreddistfin3 =
    ypreddistfin3 +
    weights_2[s] *
    h[[s]]$pred_dist
  
  # Weighted pointwise log-likelihoods
  elpd1aux =
    elpd1aux +
    weights[s] *
    h[[s]]$log_lik
  
  elpd3aux =
    elpd3aux +
    weights_2[s] *
    h[[s]]$log_lik
  
  # elpdind[s] =
  #   mean(
  #     apply(
  #       h[[s]]$log_lik,
  #       2,
  #       sum
  #     )
  #   )
}


############################################################
# ELPD for benchmark models and model averages
############################################################

elpdsort =
  sort(
    elpdind,
    decreasing = TRUE
  )


# PSIS-LOO for the model based on the reported spatial graph
elpdrepaux =
  loo(
    res_rep$log_lik
  )

elpdrep =
  elpdrepaux$estimates[
    "elpd_loo",
    "Estimate"
  ]

seelpdrep =
  elpdrepaux$estimates[
    "elpd_loo",
    "SE"
  ]


# PSIS-LOO for the temporal-only model
elpdtimeaux =
  loo(
    res_time$log_lik
  )

elpdtime =
  elpdtimeaux$estimates[
    "elpd_loo",
    "Estimate"
  ]

seelpdtime =
  elpdtimeaux$estimates[
    "elpd_loo",
    "SE"
  ]


# PSIS-LOO for the randomly selected graph
elpdrandaux =
  loo(
    res_rand$log_lik
  )

elpdrand =
  elpdrandaux$estimates[
    "elpd_loo",
    "Estimate"
  ]

seelpdrand =
  elpdrandaux$estimates[
    "elpd_loo",
    "SE"
  ]


# PSIS-LOO calculated from the stacking-weighted log-likelihood
elpd1aux2 =
  loo(
    elpd1aux
  )

elpd1 =
  elpd1aux2$estimates[
    "elpd_loo",
    "Estimate"
  ]

seelpd1 =
  elpd1aux2$estimates[
    "elpd_loo",
    "SE"
  ]


# PSIS-LOO calculated from the pseudo-BMA weighted log-likelihood
elpd3aux2 =
  loo(
    elpd3aux
  )

elpd3 =
  elpd3aux2$estimates[
    "elpd_loo",
    "Estimate"
  ]

seelpd3 =
  elpd3aux2$estimates[
    "elpd_loo",
    "SE"
  ]


# Best individual candidate ELPD
elpdsort[1]


############################################################
# Display predictive performance summaries
############################################################

#### Reported spatial graph
elpdrep
seelpdrep

#### Stacking
elpd1
seelpd1

#### Pseudo-BMA
elpd3
seelpd3

#### Best candidate tree
elpdindord[1]
seord[1]

#### Random tree
elpdrand
seelpdrand

#### No spatial dependence
elpdtime
seelpdtime


############################################################
# Identify the best candidate spatial graph
############################################################

# Candidate with the largest individual ELPD.
bestmod = h[[which(elpdind == max(elpdsort))]]

# Posterior predictive draws from the best candidate model
bestpred =
  bestmod$pred_dist


############################################################
# Out-of-sample predictive evaluation
############################################################

# True observations for the forecasting period
yreal =
  datapred$logTMA


############################################################
# Temporal-only model
############################################################

# Posterior 2.5%, 50%, and 97.5% predictive quantiles
ictime =
  t(
    apply(
      res_time$pred_dist,
      2,
      quantile,
      probs = c(
        0.025,
        0.5,
        0.975
      )
    )
  )

# Coverage indicator for each future observation
inditime =
  as.numeric(
    ictime[, 1] <= yreal &
      ictime[, 3] >= yreal
  )

# Empirical coverage probability
porctime =
  sum(inditime) /
  length(inditime)

# RMSE based on posterior median predictions
rmsetime =
  sqrt(
    mean(
      (
        ictime[, 2] -
          datapred$logTMA
      )^2
    )
  )

# Average width of the predictive intervals
lengthtime =
  mean(
    ictime[, 3] -
      ictime[, 1]
  )


############################################################
# Random spatial graph
############################################################

icrand =
  t(
    apply(
      res_rand$pred_dist,
      2,
      quantile,
      probs = c(
        0.025,
        0.5,
        0.975
      )
    )
  )

indirand =
  as.numeric(
    icrand[, 1] <= yreal &
      icrand[, 3] >= yreal
  )

porcrand =
  sum(indirand) /
  length(indirand)

rmserand =
  sqrt(
    mean(
      (
        icrand[, 2] -
          datapred$logTMA
      )^2
    )
  )

lengthrand =
  mean(
    icrand[, 3] -
      icrand[, 1]
  )


############################################################
# Reported spatial graph
############################################################

icrep =
  t(
    apply(
      res_rep$pred_dist,
      2,
      quantile,
      probs = c(
        0.025,
        0.5,
        0.975
      )
    )
  )

indirep =
  as.numeric(
    icrep[, 1] <= yreal &
      icrep[, 3] >= yreal
  )

porcrep =
  sum(indirep) /
  length(indirep)

rmserep =
  sqrt(
    mean(
      (
        icrep[, 2] -
          datapred$logTMA
      )^2
    )
  )

lengthrep =
  mean(
    icrep[, 3] -
      icrep[, 1]
  )


############################################################
# Stacking predictions
############################################################

icstack =
  t(
    apply(
      ypreddistfin,
      2,
      quantile,
      probs = c(
        0.025,
        0.5,
        0.975
      )
    )
  )

indistack =
  as.numeric(
    icstack[, 1] <= yreal &
      icstack[, 3] >= yreal
  )

porcstack =
  sum(indistack) /
  length(indistack)

rmsestack =
  sqrt(
    mean(
      (
        icstack[, 2] -
          datapred$logTMA
      )^2
    )
  )

lengthstack =
  mean(
    icstack[, 3] -
      icstack[, 1]
  )


############################################################
# Pseudo-BMA predictions
############################################################

icbma =
  t(
    apply(
      ypreddistfin3,
      2,
      quantile,
      probs = c(
        0.025,
        0.5,
        0.975
      )
    )
  )

indibma =
  as.numeric(
    icbma[, 1] <= yreal &
      icbma[, 3] >= yreal
  )

porcbma =
  sum(indibma) /
  length(indibma)

rmsebma =
  sqrt(
    mean(
      (
        icbma[, 2] -
          datapred$logTMA
      )^2
    )
  )

lengthbma =
  mean(
    icbma[, 3] -
      icbma[, 1]
  )


############################################################
# Best candidate graph
############################################################

icbest =
  t(
    apply(
      bestpred,
      2,
      quantile,
      probs = c(
        0.025,
        0.5,
        0.975
      )
    )
  )

indibest =
  as.numeric(
    icbest[, 1] <= yreal &
      icbest[, 3] >= yreal
  )

porcbest =
  sum(indibest) /
  length(indibest)

rmsebest =
  sqrt(
    mean(
      (
        icbest[, 2] -
          datapred$logTMA
      )^2
    )
  )

lengthbest =
  mean(
    icbest[, 3] -
      icbest[, 1]
  )


############################################################
# First candidate graph / WAIC-selected graph
############################################################

# This evaluates the first element of h. In the original script,
# this appears to correspond to the graph selected using WAIC.
icbestwaic =
  t(
    apply(
      h[[1]]$pred_dist,
      2,
      quantile,
      probs = c(
        0.025,
        0.5,
        0.975
      )
    )
  )

indibestwaic =
  as.numeric(
    icbestwaic[, 1] <= yreal &
      icbestwaic[, 3] >= yreal
  )

porcbestwaic =
  sum(indibestwaic) /
  length(indibestwaic)

rmsebestwaic =
  sqrt(
    mean(
      (
        icbestwaic[, 2] -
          datapred$logTMA
      )^2
    )
  )

lengthbestwaic =
  mean(
    icbestwaic[, 3] -
      icbestwaic[, 1]
  )


############################################################
# Collect performance measures
############################################################

# ELPD values
elpds =
  c(
    elpdrep,
    elpdrand,
    elpdsort[1],
    elpd1,
    elpd3,
    elpdtime
  )

# Root mean squared prediction errors
rmses =
  c(
    rmserep,
    rmserand,
    rmsebest,
    rmsestack,
    rmsebma,
    rmsetime
  )

# Empirical predictive interval coverage
porcs =
  c(
    porcrep,
    porcrand,
    porcbest,
    porcstack,
    porcbma,
    porctime
  )

# Mean predictive interval widths
lengths =
  c(
    lengthrep,
    lengthrand,
    lengthbest,
    lengthstack,
    lengthbma,
    lengthtime
  )



#########################
# Criteria Summary table#
#########################

summ = cbind(
  elpds,
  rmses,
  porcs,
  lengths
)

rownames(summ) = c("Reported","Random tree", "Best", "Stacking","BMA","No_spatial")

summ

############################################################
# Pairwise PSIS-LOO comparisons
############################################################

# PSIS-LOO object for the best candidate graph
elpdbest =
  loo(
    bestmod$log_lik
  )

# Compare the best candidate against each alternative model
loo_compare(
  elpdbest,
  elpd1aux2
)

loo_compare(
  elpdbest,
  elpd3aux2
)

loo_compare(
  elpdbest,
  elpdrandaux
)

loo_compare(
  elpdbest,
  elpdtimeaux
)

loo_compare(
  elpdbest,
  elpdrepaux
)


############################################################
# Prepare data for spatial visualization
############################################################

library(ggplot2)
library(dplyr)

# Ensure observations are ordered consistently by zone and time
datapred <-
  datapred %>%
  arrange(
    zone,
    t
  )


############################################################
# Spatial maps
############################################################

library(sf)
library(dplyr)
library(ggplot2)


# Read shapefile containing the spatial units
chile_shp <-
  st_read(
    "chile_shp_union.shp"
  )


# Convert zone to character to facilitate joining with the shapefile,
# and ensure that the date variable is stored as Date.
datapred <-
  datapred %>%
  mutate(
    zone =
      as.character(zone),
    
    date =
      as.Date(date)
  )



############################################################
# Aggregate predictions by spatial zone and forecasting date
############################################################

pred_mapa <-
  datapred %>%
  group_by(
    date,
    zone
  ) %>%
  summarise(
    
    # Best candidate graph prediction
    ypred_map =
      mean(
        ypredbest,
        na.rm = TRUE
      ),
    
    # Stacking prediction
    ypredstack_map =
      mean(
        ypred,
        na.rm = TRUE
      ),
    
    # Upper predictive limit for the best graph
    ypred975_map =
      mean(
        ypredbest975,
        na.rm = TRUE
      ),
    
    # Observed log response
    logTMA_map =
      mean(
        logTMA,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


############################################################
# Join spatial polygons with prediction summaries
############################################################

mapa_pred <-
  chile_shp %>%
  left_join(
    pred_mapa,
    by = "zone"
  )


############################################################
# Faceted prediction map: best candidate graph
############################################################

# Display predictions for all seven forecasting dates
# using one map per date.
prediction <-
  ggplot(mapa_pred) +
  
  geom_sf(
    aes(
      fill = ypred_map
    ),
    color = "black",
    linewidth = 0.2
  ) +
  
  facet_wrap(
    ~ date,
    ncol = 2,
    nrow = 4
  ) +
  
  labs(
    fill = "",
    title = ""
  ) +
  
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0
  ) +
  
  coord_sf(
    expand = FALSE,
    default_crs =
      sf::st_crs(4326)
  ) +
  
  scale_x_continuous(
    breaks =
      c(
        -73.14,
        -73.10,
        -73.06
      ),
    
    labels =
      function(x)
        paste0(
          abs(x),
          "°W"
        )
  ) +
  
  scale_y_continuous(
    breaks =
      c(
        -36.89,
        -36.86,
        -36.83
      ),
    
    labels =
      function(y)
        paste0(
          abs(y),
          "°S"
        )
  ) +
  
  theme_minimal()


# Range of the predicted values
lims1 <-
  range(
    mapa_pred$ypred_map,
    na.rm = TRUE
  )


# Forecasting dates represented in the maps
fechas =
  unique(
    mapa_pred$date
  )


############################################################
# Faceted map of observed values
############################################################

# Observed values from 2021-08-11 to 2021-08-17
real <-
  ggplot(mapa_pred) +
  
  geom_sf(
    aes(
      fill = logTMA_map
    ),
    color = "black",
    linewidth = 0.2
  ) +
  
  facet_wrap(
    ~ date,
    ncol = 2,
    nrow = 4
  ) +
  
  labs(
    fill = "",
    title = ""
  ) +
  
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0
  ) +
  
  coord_sf(
    expand = FALSE,
    default_crs =
      sf::st_crs(4326)
  ) +
  
  scale_x_continuous(
    breaks =
      c(
        -73.14,
        -73.10,
        -73.06
      ),
    
    labels =
      function(x)
        paste0(
          abs(x),
          "°W"
        )
  ) +
  
  scale_y_continuous(
    breaks =
      c(
        -36.89,
        -36.86,
        -36.83
      ),
    
    labels =
      function(y)
        paste0(
          abs(y),
          "°S"
        )
  ) +
  
  theme_minimal()


############################################################
# Spatial prediction errors: best candidate graph
############################################################

# Range of prediction errors for the best candidate graph
lims2 <-
  range(
    mapa_pred$ypred_map -
      mapa_pred$logTMA_map,
    na.rm = TRUE
  )


# Difference between the best-graph prediction
# and the observed value.
diff <-
  ggplot(mapa_pred) +
  
  geom_sf(
    aes(
      fill =
        ypred_map -
        logTMA_map
    ),
    color = "black",
    linewidth = 0.2
  ) +
  
  facet_wrap(
    ~ date,
    ncol = 2,
    nrow = 4
  ) +
  
  labs(
    fill = "",
    title = ""
  ) +
  
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0
  ) +
  
  coord_sf(
    expand = FALSE,
    default_crs =
      sf::st_crs(4326)
  ) +
  
  scale_x_continuous(
    breaks =
      c(
        -73.14,
        -73.10,
        -73.06
      ),
    
    labels =
      function(x)
        paste0(
          abs(x),
          "°W"
        )
  ) +
  
  scale_y_continuous(
    breaks =
      c(
        -36.89,
        -36.86,
        -36.83
      ),
    
    labels =
      function(y)
        paste0(
          abs(y),
          "°S"
        )
  ) +
  
  theme_minimal()


############################################################
# Spatial prediction errors: stacking
############################################################

diff_stack <-
  ggplot(mapa_pred) +
  
  geom_sf(
    aes(
      fill =
        ypredstack_map -
        logTMA_map
    ),
    color = "black",
    linewidth = 0.2
  ) +
  
  facet_wrap(
    ~ date,
    ncol = 2,
    nrow = 4
  ) +
  
  labs(
    fill = "",
    title = ""
  ) +
  
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    
    # Use the same error scale as the best-graph map
    # to facilitate direct comparison.
    limits = lims2
  ) +
  
  coord_sf(
    expand = FALSE,
    default_crs =
      sf::st_crs(4326)
  ) +
  
  scale_x_continuous(
    breaks =
      c(
        -73.14,
        -73.10,
        -73.06
      ),
    
    labels =
      function(x)
        paste0(
          abs(x),
          "°W"
        )
  ) +
  
  scale_y_continuous(
    breaks =
      c(
        -36.89,
        -36.86,
        -36.83
      ),
    
    labels =
      function(y)
        paste0(
          abs(y),
          "°S"
        )
  ) +
  
  theme_minimal()


# Display the best-graph prediction map
prediction


############################################################
# Save faceted figures
############################################################

ggsave(
  "diff_stack.eps",
  plot = diff_stack,
  device = cairo_ps,
  dpi = 600,
  width = 8,
  height = 6,
  units = "in"
)


ggsave(
  "diffbest.eps",
  plot = diff,
  device = cairo_ps,
  dpi = 600,
  width = 8,
  height = 6,
  units = "in"
)


# NOTE:
# In the original script, prediction_stack is saved here,
# but it has not yet been defined at this point.
# It is created later inside the date-specific loop.
ggsave(
  "prediction_stack.eps",
  plot = prediction_stack,
  device = cairo_ps,
  dpi = 600,
  width = 8,
  height = 6,
  units = "in"
)


ggsave(
  "prediction_best.eps",
  plot = prediction,
  device = cairo_ps,
  dpi = 600,
  width = 8,
  height = 6,
  units = "in"
)


ggsave(
  "realdata.eps",
  plot = real,
  device = cairo_ps,
  dpi = 600,
  width = 8,
  height = 6,
  units = "in"
)


############################################################
# Generate one stacking prediction map per forecasting date
############################################################

# Observed range, retained for reference
lims3 =
  range(
    mapa_pred$logTMA_map,
    na.rm = TRUE
  )


for (i in 1:length(fechas)) {
  
  # Stacking prediction for forecasting date i
  prediction_stack <-
    ggplot(
      mapa_pred %>%
        filter(
          date == fechas[i]
        )
    ) +
    
    labs(
      fill = "",
      title = fechas[i]
    ) +
    
    geom_sf(
      aes(
        fill =
          ypredstack_map
      ),
      color = "black",
      linewidth = 0.2
    ) +
    
    # A common scale is used across dates so that
    # spatial and temporal patterns are directly comparable.
    scale_fill_gradientn(
      colours =
        c(
          "#FEE5D9",
          "#FCAE91",
          "#FB6A4A",
          "#DE2D26",
          "#A50F15"
        ),
      
      limits =
        c(
          1.5,
          3.2
        ),
      
      breaks =
        c(
          2,
          2.5,
          3
        ),
      
      oob =
        scales::squish
    ) +
    
    coord_sf(
      expand = FALSE,
      default_crs =
        sf::st_crs(4326)
    ) +
    
    scale_x_continuous(
      breaks =
        c(
          -73.14,
          -73.10,
          -73.06
        ),
      
      labels =
        function(x)
          paste0(
            abs(x),
            "°W"
          )
    ) +
    
    scale_y_continuous(
      breaks =
        c(
          -36.89,
          -36.86,
          -36.83
        ),
      
      labels =
        function(y)
          paste0(
            abs(y),
            "°S"
          )
    ) +
    
    theme_minimal() +
    
    theme(
      plot.title =
        element_text(
          size = 22,
          face = "bold"
        ),
      
      axis.title =
        element_text(
          size = 16
        ),
      
      axis.text =
        element_text(
          size = 14
        ),
      
      legend.title =
        element_text(
          size = 16
        ),
      
      legend.text =
        element_text(
          size = 14
        )
    )
  
  
  # Directory where individual maps are stored
  out_dir =
    "C:/Users/alejo/Dropbox/postdoc20242/random_graphs_new/Application COVID/mapas/"
  
  
  ggsave(
    paste0(
      out_dir,
      "pred_stack_",
      i,
      ".eps"
    ),
    plot = prediction_stack,
    device = cairo_ps,
    dpi = 600,
    width = 8,
    height = 6,
    units = "in"
  )
}


############################################################
# Generate one best-graph prediction map per forecasting date
############################################################

for (i in 1:length(fechas)) {
  
  prediction <-
    ggplot(
      mapa_pred %>%
        filter(
          date == fechas[i]
        )
    ) +
    
    labs(
      fill = "",
      title = fechas[i]
    ) +
    
    geom_sf(
      aes(
        fill =
          ypred_map
      ),
      color = "black",
      linewidth = 0.2
    ) +
    
    scale_fill_gradientn(
      colours =
        c(
          "#FEE5D9",
          "#FCAE91",
          "#FB6A4A",
          "#DE2D26",
          "#A50F15"
        ),
      
      limits =
        c(
          1.5,
          3.2
        ),
      
      breaks =
        c(
          2,
          2.5,
          3
        ),
      
      oob =
        scales::squish
    ) +
    
    coord_sf(
      expand = FALSE,
      default_crs =
        sf::st_crs(4326)
    ) +
    
    scale_x_continuous(
      breaks =
        c(
          -73.14,
          -73.10,
          -73.06
        ),
      
      labels =
        function(x)
          paste0(
            abs(x),
            "°W"
          )
    ) +
    
    scale_y_continuous(
      breaks =
        c(
          -36.89,
          -36.86,
          -36.83
        ),
      
      labels =
        function(y)
          paste0(
            abs(y),
            "°S"
          )
    ) +
    
    theme_minimal() +
    
    theme(
      plot.title =
        element_text(
          size = 22,
          face = "bold"
        ),
      
      axis.title =
        element_text(
          size = 16
        ),
      
      axis.text =
        element_text(
          size = 14
        ),
      
      legend.title =
        element_text(
          size = 16
        ),
      
      legend.text =
        element_text(
          size = 14
        )
    )
  
  
  out_dir =
    "C:/Users/alejo/Dropbox/postdoc20242/random_graphs_new/Application COVID/mapas/"
  
  
  ggsave(
    paste0(
      out_dir,
      "pred_best_",
      i,
      ".eps"
    ),
    plot = prediction,
    device = cairo_ps,
    dpi = 600,
    width = 8,
    height = 6,
    units = "in"
  )
}


############################################################
# Generate one observed-value map per forecasting date
############################################################

for (i in 1:length(fechas)) {
  
  real <-
    ggplot(
      mapa_pred %>%
        filter(
          date == fechas[i]
        )
    ) +
    
    labs(
      fill = "",
      title = fechas[i]
    ) +
    
    geom_sf(
      aes(
        fill =
          logTMA_map
      ),
      color = "black",
      linewidth = 0.2
    ) +
    
    # Use the same scale as the prediction maps
    # to make visual comparisons meaningful.
    scale_fill_gradientn(
      colours =
        c(
          "#FEE5D9",
          "#FCAE91",
          "#FB6A4A",
          "#DE2D26",
          "#A50F15"
        ),
      
      limits =
        c(
          1.5,
          3.2
        ),
      
      breaks =
        c(
          2,
          2.5,
          3
        ),
      
      oob =
        scales::squish
    ) +
    
    coord_sf(
      expand = FALSE,
      default_crs =
        sf::st_crs(4326)
    ) +
    
    scale_x_continuous(
      breaks =
        c(
          -73.14,
          -73.10,
          -73.06
        ),
      
      labels =
        function(x)
          paste0(
            abs(x),
            "°W"
          )
    ) +
    
    scale_y_continuous(
      breaks =
        c(
          -36.89,
          -36.86,
          -36.83
        ),
      
      labels =
        function(y)
          paste0(
            abs(y),
            "°S"
          )
    ) +
    
    theme_minimal() +
    
    theme(
      plot.title =
        element_text(
          size = 22,
          face = "bold"
        ),
      
      axis.title =
        element_text(
          size = 16
        ),
      
      axis.text =
        element_text(
          size = 14
        ),
      
      legend.title =
        element_text(
          size = 16
        ),
      
      legend.text =
        element_text(
          size = 14
        )
    )
  
  
  out_dir =
    "C:/Users/alejo/Dropbox/postdoc20242/random_graphs_new/Application COVID/mapas/"
  
  
  ggsave(
    paste0(
      out_dir,
      "realdata_",
      i,
      ".eps"
    ),
    plot = real,
    device = cairo_ps,
    dpi = 600,
    width = 8,
    height = 6,
    units = "in"
  )
}


############################################################
# Generate stacking prediction-error maps by date
############################################################

for (i in 1:length(fechas)) {
  
  # Maximum absolute stacking prediction error.
  # This quantity is computed here but is not subsequently
  # used in the current plotting scale.
  lim_diff <-
    max(
      abs(
        mapa_pred$logTMA_map -
          mapa_pred$ypredstack_map
      ),
      na.rm = TRUE
    )
  
  
  # Error is defined here as observed value minus prediction.
  diff_stack <-
    ggplot(
      mapa_pred %>%
        filter(
          date == fechas[i]
        )
    ) +
    
    labs(
      fill = "",
      title = fechas[i]
    ) +
    
    geom_sf(
      aes(
        fill =
          logTMA_map -
          ypredstack_map
      ),
      color = "black",
      linewidth = 0.2
    ) +
    
    scale_fill_gradientn(
      colours =
        c(
          "#FEE5D9",
          "#FCAE91",
          "#FB6A4A",
          "#DE2D26",
          "#A50F15"
        ),
      
      limits =
        c(
          -1.7,
          0.8
        ),
      
      breaks =
        c(
          -1,
          -0.5,
          0
        ),
      
      oob =
        scales::squish
    ) +
    
    coord_sf(
      expand = FALSE,
      default_crs =
        sf::st_crs(4326)
    ) +
    
    scale_x_continuous(
      breaks =
        c(
          -73.14,
          -73.10,
          -73.06
        ),
      
      labels =
        function(x)
          paste0(
            abs(x),
            "°W"
          )
    ) +
    
    scale_y_continuous(
      breaks =
        c(
          -36.89,
          -36.86,
          -36.83
        ),
      
      labels =
        function(y)
          paste0(
            abs(y),
            "°S"
          )
    ) +
    
    theme_minimal() +
    
    theme(
      plot.title =
        element_text(
          size = 22,
          face = "bold"
        ),
      
      axis.title =
        element_text(
          size = 16
        ),
      
      axis.text =
        element_text(
          size = 14
        ),
      
      legend.title =
        element_text(
          size = 16
        ),
      
      legend.text =
        element_text(
          size = 14
        )
    )
  
  
  out_dir =
    "C:/Users/alejo/Dropbox/postdoc20242/random_graphs_new/Application COVID/mapas/"
  
  
  ggsave(
    paste0(
      out_dir,
      "diff_stack_",
      i,
      ".eps"
    ),
    plot = diff_stack,
    device = cairo_ps,
    dpi = 600,
    width = 8,
    height = 6,
    units = "in"
  )
}


############################################################
# Generate best-graph prediction-error maps by date
############################################################

for (i in 1:length(fechas)) {
  
  # Error is defined as observed value minus
  # the prediction from the best candidate graph.
  diff <-
    ggplot(
      mapa_pred %>%
        filter(
          date == fechas[i]
        )
    ) +
    
    labs(
      fill = "",
      title = fechas[i]
    ) +
    
    geom_sf(
      aes(
        fill =
          logTMA_map -
          ypred_map
      ),
      color = "black",
      linewidth = 0.2
    ) +
    
    # The same limits are used as in the stacking error maps
    # to allow direct comparison between the two approaches.
    scale_fill_gradientn(
      colours =
        c(
          "#FEE5D9",
          "#FCAE91",
          "#FB6A4A",
          "#DE2D26",
          "#A50F15"
        ),
      
      limits =
        c(
          -1.7,
          0.8
        ),
      
      breaks =
        c(
          -1,
          -0.5,
          0
        ),
      
      oob =
        scales::squish
    ) +
    
    coord_sf(
      expand = FALSE,
      default_crs =
        sf::st_crs(4326)
    ) +
    
    scale_x_continuous(
      breaks =
        c(
          -73.14,
          -73.10,
          -73.06
        ),
      
      labels =
        function(x)
          paste0(
            abs(x),
            "°W"
          )
    ) +
    
    scale_y_continuous(
      breaks =
        c(
          -36.89,
          -36.86,
          -36.83
        ),
      
      labels =
        function(y)
          paste0(
            abs(y),
            "°S"
          )
    ) +
    
    theme_minimal() +
    
    theme(
      plot.title =
        element_text(
          size = 22,
          face = "bold"
        ),
      
      axis.title =
        element_text(
          size = 16
        ),
      
      axis.text =
        element_text(
          size = 14
        ),
      
      legend.title =
        element_text(
          size = 16
        ),
      
      legend.text =
        element_text(
          size = 14
        )
    )
  
  
  out_dir =
    "C:/Users/alejo/Dropbox/postdoc20242/random_graphs_new/Application COVID/mapas/"
  
  
  ggsave(
    paste0(
      out_dir,
      "diff_best_",
      i,
      ".eps"
    ),
    plot = diff,
    device = cairo_ps,
    dpi = 600,
    width = 8,
    height = 6,
    units = "in"
  )
}