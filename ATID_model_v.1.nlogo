;********************************************************************;
;                      GENERAL INFORMATION                           ;
;********************************************************************;
;This is a trait- and individual-based model to investigate,
; which environmental factors and plant functional traits are important for plant performance and community dynamics in deserts.
; The model is an integral part of the paper Zakharova L., Meyer K.M., Seifan M. (20XX) Combining trait- and individual-based modelling to understand desert plant community dynamics.
;
; SLIDERS
; in-s-number-an & in-s-number-ma control the number of seeds of each species initialized in the model.
; in-ad-number-an & in-ad-number-ma control the number of adult plants of each species initialized in the model.
; The postfixes "-an" & "-ma" stand for two generic species with contrasting dispersal strategies, after Anastatica hierochuntica & Malva parviflora.
; TOPOGRAPHY group:
; min-height defines the height of the lowest patch, controls the steepness of the site
; height-dif defines the difference in height between the lowest patch and the highest patch, controls the steepness of the site
; sd-height defines standard deviation from the mean of the absolutely smooth slope, controls the smoothness of the site
; PRECIPITATION group:
; min-days-year defines the minimal number of rainy days in the season possible for the specific site
; rain-days-year-dif defines the difference between the minimal and the maximal number of rainy days in the season possible for the specific site
; rain-season defines typical length of the vegetation season for the specific site
; dew-days-before defines the number of days before the first rain when dew deposition might affect plant growth
; dew-days-after defines the number of days after the last rain when dew deposition might affect plant growth
; WATER AVAILABLE FOR PLANTS group
; av-water-min controls how much water is available for plants at the topographical object “local minimum” after a rain event
; av-water-dif-min-slope defines how much less water is available for plants at the topographical object “slope” after a rain event in comparison to the topographical object “local minimum”
; av-water-dif-slope-max defines how much less water is available for plants at the topographical object “local maximum” after a rain event in comparison to the topographical object “slope”
; av-water-dif-max-dew defines how much less water is available for plants at all topographical objects as a result of dew deposition
; in comparison to the topographical object “local maximum”
; Note: all parameters defining the difference between the parameters are introduced for purposes of sensitivity analysis.
; For the model version for simulations, they should be replaced with the sliders controlling the parameter values directly.
; Anastatica & Malva traits control sets of plant functional traits for two species.
;
; SWITCHES
; dew? chooses between the scenario with dew and without it.
; zoi-approach chooses the algorithm calculating zones-of-influence and plant growth following either Lin et al. 2012 or Weiner & Damgaard, 2006.
; References:
; 1. Lin, Y., Berger, U., Grimm, V., & Ji, Q. R. (2012). Differences between symmetric and asymmetric facilitation matter: exploring the interplay between modes of
; positive and negative plant interactions. Journal of Ecology, 100(6), 1482-1491.
; 2. Weiner, J., & Damgaard, C. (2006). Size-asymmetric competition and size-asymmetric growth in a spatially explicit zone-of-influence model of plant competition.
; Ecological Research, 21(5), 707-712.


;********************************************************************;
;                     Agent types, actors                            ;
;********************************************************************;

globals
[
  size-factor       ; Takes care of size depending on cell size, after Lin et al., 2012
  days              ; Number of days, equal to ticks
  years             ; Number of years
  ; Storage variables, used for output
  ; Current amount of adults of every species
  ad-number-An      ; Number of adult Anastatica
  ad-number-Ma      ; Number of adult Malva
  ; Current amount of seeds of every species
  s-number-An       ; Number of seeds, Anastatica
  s-number-Ma       ; Number of seeds, Malva
  bio-mean-an       ; Mean biomass of all adult Anastatica
  bio-mean-ma       ; Mean biomass of all adult Malva
  rainy-days        ; A list with the numbers of days with a rain event during the year
  first-rain-day    ; The first rainy day in a season, defined as the first number in the list of simulated rainy days
  last-rain-day     ; The last rainy day in a season, defined as the last number in the list of simulated rainy days
  season-correction ; Number of days added at the beginning of the vegetation season to simulate the effect of dew on its length
]

;********************************************************************;
; AGENTS' AND PATCHES' CREATION:                                     ;
; Patches: topography & water availability (dependent on topography) ;
; Agents: two species with a set of traits (adult plants + seeds)    ;
;********************************************************************;
; plants
breed [plants plant]

; patches are 1x1 m, soil salinity affects soil water availability and isn't described explicitly
; soil water availability depends also on topography and soil properties,
; in the current version it is described as a probability that plants have access to soil water on the specific patch
patches-own
[
p-biomass_acc          ; Accumulated biomass of competing plants, used in the procedure Competition
p-biomass_acc_inv      ; Inverted value of accumulated biomass
p-water-availability   ; soil water availability, a resulting parameter of precipitation of all kinds, water content, salinity and topography; the probability that plants obtain soil water
height                 ; Height above the sea level; the user can define the range of the lowest (left, down patch -100, -50) and
                       ; the highest points as well as sd from a mean for each patch;
                       ; in the version for sensitivity analysis, the parameters are the lowest height,
                       ; the difference between the lowest and the highest and standard deviation from a mean, which is defined as a slope
act-height             ; This is a variable recording the height of the patch of interest to adjust heights of its neighbours accordingly; part of topography simulation
height-fixed?          ; This variable prevents setting patches' height more than once; part of topography simulation
topography-fixed?      ; This variable prevents setting patches' topography more than once; part of topography simulation
topography             ; The topography defined based on simulated heights above the sea level
                       ; "slope" is a patch which has both neighbouring patches lying above and below it
                       ; "local maximum" is a patch surrounded by the patches lying below it
                       ; "local minimum" is a patch surrounded by the patches lying above it
                       ; depending on the topography there are different rules for catching dew, collecting rainwater
                       ; and consequently, different probabilities for plants to grow and for seeds to germinate
]

;********************************************************************;
;          Individual parameters and variables of the agents         ;
;********************************************************************;
plants-own
[
 adult?                   ; Define if there is a seed or an adult plant
 profile                  ; Species identity with a set of traits
 ad-biomass               ; Current biomass, for seeds, equals s-mass
 ad-biomass-max           ; Maximum biomass a species can achieve, equal to the asymptotic body size
 age                      ; Time since germination for plants and since creation for seeds, in days
 ad-prob-surv             ; Probability of a plant to survive as an entire organism
 juv-biomass              ; A fraction of maximum biomass which determines the threshold between juvenile and adult plants
 juv-surv                 ; Coefficient reducing the juvenile probability to survive relative to the adult probability
  ; Growth
 ad-rgr                   ; Relative growth rate, part of the current biomass equations after Lin et al. 2012 or Weiner & Damgaard, 2006.
  ; Reproduction
 ad-age-repr              ; Age at that a plant can start reproducing under suitable conditions
 ad-biomass-repr          ; Reproduction biomass, a threshold which plants should reach to start the production of seeds
 s-mass                   ; Seed weight, which lies in a species-specific range
 s-sp-number              ; Species-specific number of seeds throughout the entire vegetation period
 s-number-actual          ; Number of seeds to be produced if a plant reproduces depending on current biomass
 count-s-number           ; Counter of produced seeds, prevents overproduction of seeds
 ; Dispersal
 s-disp-type              ; Defines to which dispersal type a plant belongs to:
                          ; escape (seeds are produced and almost immediately distributed) or
                          ; protection (seeds are preserved until the conditions are suitable for dispersion/ germination) (Guttermann, 2000)
 s-disp-age               ; Age at what a seed is dispersed after the creation
 s-dispersed?             ; This variable shows status if a seed is still on the parent plant
 s-disp-scale             ; Species-specific parameter δ in Weibull-shaped dispersal; defines how far from a parent plant seeds will be brought
 s-disp-shape             ; Species-specific parameter β in Weibull-shaped dispersal
  ; Germination
 s-days-dorm              ; Number of days that seeds should spend in dormancy; defines the possible start of germination
 s-pr-germ                ; Probability of seed germination, if other conditions are met
 ; Competition
 ad-radius                ; The radius of zone-of-influence (ZOI), after Lin et al. 2012
 ad-zoi-patches           ; ZOI: number of patches inside the radius of ZOI
 ad-zoi-overlap           ; ZOI: effective area, calculated for every plant as a difference between the area the plant covers and the area lost to competition with its neighbours
 ad-comp-asymmetry        ; Degree of asymmetry, if it equals 0, resources are shared equally, regardless of their biomass.
                          ; if it equals 1, resources are distributed proportionally to the biomass (perfect size-symmetry) (Weiner & Damgaard, 2006)
]



;********************************************************************;
;                           INITIALIZATION                           ;
;********************************************************************;

to setup-environment
  clear-all               ; Clean-up the grid
  reset-ticks             ; and the time from former run
  set years 1
  set days 1
  set-default-shape plants "circle 2"
  set size-factor  ( world-width / 5000) ; The world is 100*50 m, each patch - 1 m^2.

;********************************************************************;
;                       Topography simulation                        ;
;********************************************************************;

; In case of real-world versions of the model, this part of the model is replaced by values from a digital elevation model.
; Set the general slope and the height of the site
; The slope is created by setting up one corner as the lowest site and by gradually simulating the height distribution up to the highest one
; In the version of the model for sensitivity analysis, the largest height is substituted by (min-height + height-dif).

ask patches
  [
    set height-fixed? FALSE
  ]
; The lowest height is assigned to the left lower corner
ask patch min-pxcor min-pycor
  [
    set height min-height
    set act-height height
; After the height of the patch is changed, the patch is excluded from the following changes
    set height-fixed? TRUE
  ]
ask patches
  [
    ; Those patches are addressed that have already their height set up
   ask patches with [height = act-height AND height-fixed? = TRUE]
    [
     let my-pxcor pxcor
     let my-pycor pycor
     let my-height height
      ; Those patches are addressed that are lying to the right and above the initial lowest patch
     ask neighbors with [pxcor <= my-pxcor + 1  AND pycor <= my-pycor + 1 AND height-fixed? = FALSE ]
       [
        set height my-height + random-normal (height-dif / world-width) sd-height
; The increment is defined, by which the neighbouring patches can be higher than the lowest neighbour;
; this algorithm leads to a gradual increase in height, forming a slope. sd-height defines the smoothness of the slope
; The simulated site height should not extend the max-height according to the measurements in the reality
          if height > (min-height + height-dif) [set height (min-height + height-dif)]
; The same steps are repeated for the neighbouring patch as for the lowest patch so that they are excluded from the process of setting up the height
        set act-height height
        set height-fixed? TRUE
          let col-height int (height / 10)
        set pcolor col-height
       ]
    ]
 ; Now different topographic properties are assigned to the patches depending on their location relative to their neighbours
 ; All patches are slopes by default
    set topography "slope"
    set topography-fixed? FALSE
 ; To avoid asking again already set minima & maxima
    ask patches with [height-fixed? = TRUE AND topography = "slope" and topography-fixed? = FALSE]
    [
    set topography-fixed? TRUE
    let my-height height
 ; Patches are counted strictly below the patch of interest
    let number-below-neighbors count neighbors with [height < my-height]
 ; Patches are counted strictly above the patch of interest
    let number-above-neighbors count neighbors with [height > my-height]
 ; If all the neighbours have height value larger or equal than this one, then this is a local minimum
    if number-below-neighbors = 0
      [
        set topography "local minimum" ; Topographical "patches" lying below or at the same height as all other in the surrounding
      ]
 ; If there is no neighbour with the height larger than this one, then there is a local maximum (but can be with equal height)
      if number-above-neighbors = 0
      [
        set topography "local maximum"  ; Patches above or at the same height as all other in the surrounding
      ]
    ]
  ]

;********************************************************************;
;                     Creation of plants and seeds                   ;
;********************************************************************;
; Assign sets of plant (functional) traits for different species
; The same species traits are assigned to both seeds and adult plants
create-plants in-s-number-an + in-ad-number-an
  [
    set profile "Anastatica"
    set adult? false ; Seeds are created
    setxy random-xcor random-ycor
    set color yellow
    set s-mass s-mass-an
    set ad-rgr ad-rgr-an
    set ad-biomass s-mass
    set age 1
    set ad-biomass-max ad-biomass-max-an
    set ad-prob-surv ad-prob-surv-an
    set juv-biomass juv-biomass-an
    set juv-surv juv-surv-an
; Reproduction
    set s-sp-number s-sp-number-an
    set count-s-number s-sp-number-an
    set ad-biomass-repr ad-biomass-repr-an
    set ad-age-repr ad-age-repr-an
; Dispersal
    ; There are two groups of plants: the first produces a lot of seeds distributed by the wind;
    ; the second are plants with serotiny; dispersal will be rain-dependent
    set s-disp-type s-disp-type-an
    set s-disp-age s-disp-age-an
    set s-dispersed? true
    set s-disp-scale s-disp-scale-an
    set s-disp-shape s-disp-shape-an
; Germination
    set s-days-dorm s-days-dorm-an
    set s-pr-germ s-pr-germ-an
; Competition
    set ad-comp-asymmetry ad-comp-asymmetry-an
    if-else zoi-approach
    [
      ; The radius of ZOI is calculated after Weiner & Damgaard, 2006
      set ad-radius sqrt ((ad-biomass ^ ( 2 / 3)) / pi)
    ]
    [
      ; The radius of ZOI is calculated after Lin et al., 2012
      set ad-radius ( ad-biomass ^ ( 3 / 8 ) ) * ( pi ^ ( -1 / 2 ) )
    ]
    ; Based on the radius, we calculate plants' apparent size
    set size ad-radius * size-factor
  ]

ask n-of in-ad-number-an plants with [profile = "Anastatica"]
  [
 ; Adult plants are created
    set adult? true
 ; Age is equal 1, for the sensitivity analysis the model is initiated with seeds only;
 ; similar to the real world at the beginning of rain (vegetation) season in the desert
    set age 1
    if-else zoi-approach
      [
      ; Current biomass is calculated after Weiner & Damgaard, 2006
        set ad-biomass s-mass + ad-rgr * (1 - (s-mass ^ 2 / ad-biomass-max ^ (4 / 3)))
      ]
      [
      ; Current biomass is calculated after Lin et al., 2012
        set ad-biomass s-mass + ad-rgr * (1 - ((s-mass / ad-biomass-max) ^ (1 / 4)))
      ]
  ]

create-plants in-s-number-ma + in-ad-number-ma
 [
   set profile "Malva"
   set adult? false
   setxy random-xcor random-ycor
   set color green
   set s-mass s-mass-ma
   set ad-rgr ad-rgr-ma
   set ad-biomass s-mass
   set age 1
   set ad-biomass-max ad-biomass-max-ma
   set ad-prob-surv ad-prob-surv-ma
   set juv-biomass juv-biomass-ma
   set juv-surv juv-surv-ma
; Reproduction
   set s-sp-number s-sp-number-ma
   set count-s-number s-sp-number-ma
   set ad-biomass-repr ad-biomass-repr-ma
   set ad-age-repr ad-age-repr-ma
; Dispersal
   set s-disp-type s-disp-type-ma
   set s-disp-age s-disp-age-ma
   set s-dispersed? true
   set s-disp-scale s-disp-scale-ma
   set s-disp-shape s-disp-shape-ma
; Germination
   set s-days-dorm s-days-dorm-ma
   set s-pr-germ s-pr-germ-ma
; Competition
   set ad-comp-asymmetry ad-comp-asymmetry-ma
    if-else zoi-approach
    [
      ; The radius of ZOI is calculated after Weiner & Damgaard, 2006
      set ad-radius sqrt ((ad-biomass ^ ( 2 / 3)) / pi)
    ]
    [
      ; The radius of ZOI is calculated after Lin et al., 2012
      set ad-radius ( ad-biomass ^ ( 3 / 8 ) ) * ( pi ^ ( -1 / 2 ) )
    ]
   ; Based on the radius, plants' apparent size is calculated
   set size ad-radius * size-factor
 ]
ask n-of in-ad-number-ma plants with [profile = "Malva"]
 [
    set adult? true
    set age 1
    if-else zoi-approach
      [
      ; Current biomass is calculated after Weiner & Damgaard, 2006
        set ad-biomass s-mass + ad-rgr * (1 - (s-mass ^ 2 / ad-biomass-max ^ (4 / 3)))
      ]
      [
      ; Current biomass is calculated after Lin et al., 2012
        set ad-biomass s-mass + ad-rgr * (1 - ((s-mass / ad-biomass-max) ^ (1 / 4)))
      ]
 ]


  ; Counting amount of the plants and seeds for OUTPUT
set ad-number-An count plants with [profile = "Anastatica" and adult? = true]
set ad-number-Ma count plants with [profile = "Malva" and adult? = true]
set s-number-An count plants with [profile = "Anastatica" and adult? = false]
set s-number-Ma count plants with [profile = "Malva" and adult? = false]

end

;********************************************************************;
; GO                                                                 ;
; RUN Model, first day of the vegetation (rain) season               ;
; Counter =1 day, Super-Counter = 1 year                             ;
;********************************************************************;

to go
if any? plants
  [
    rain                ; Environment: annually, the procedure is called once at the beginning of the year, which coincides with the vegetation season for the model
    water-availability  ; Environment: daily, changes the water availability in the soil based on either rain events or affinity of patches to accumulate water from dew
    dispersal           ; Plants: under certain conditions, depending on a type ; dispersal depends on rain events for the strategy where it caused by a raindrop
    ad-mortality        ; Plants: natural mortality - the probability to die as a part of the population (all not explicitly included factors)
    ageing              ; Plants: daily
    count-plants        ; Output
    plotting            ; Output
    set days days + 1
; Change of the years, important to organize seasonality in soil water availability alterations
    if days > 365
     [
        set years years + 1
        set days 1
     ]
    if dew? = FALSE
    [ set dew-days-before 0
      set dew-days-after 0
    ]
tick
  ]
end

to rain
 ; Rain events play a big role in seed dispersal; they also make sites in the soil depressions more attractive for growth
 ; at the beginning of every year a list of rainy days is generated, the range of rainy days and the length of the entire rainy period is defined by the user
 ; Rainy days are initialized at the beginning of the vegetation period
  if days = 1
  [
    ; A list of the rainy days in the current year,
    ; where the number of list items lies in a range of rainy days per year and the list item itself defines the serial number of rainy days in this season
    ; Alternatively, real data are used (or simulations reflecting real data) from the file
    set rainy-days n-of (min-days-year + random (min-days-year + rain-days-year-dif)) n-values rain-season [i -> i + 1 ]
    set first-rain-day min rainy-days
    set last-rain-day  max rainy-days
  ]

  ; When there is no rain for a long period -> vegetation season is over
  if-else dew? = TRUE and dew-days-before > first-rain-day
  [
    ; This correction for the vegetation season is introduced to avoid getting negative numbers
    ; when the margins of the vegetation season defined in the water-availability procedure
    set season-correction dew-days-before - first-rain-day
    set first-rain-day min rainy-days
    set last-rain-day  max rainy-days
    set rainy-days rainy-days
  ]
  [set season-correction 0]
end

; Water availability is set up as a probability that plants get water for growth
; water availability depends on the topography of the site, the difference in water availability between different topographic objects
; depends implicitly on soil properties and salinity
; Water availability is updated after rain events. It can be updated because of dew for some topographical objects if dew is considered
to water-availability
; The model runs daily (soil water availability, all main functions of plants) within the vegetation period
; Vegetation period starts with the first rain. If dew matters, then, respectively, several days earlier and finishes several days later.
  if-else days >= first-rain-day - dew-days-before + season-correction AND days <= last-rain-day + dew-days-after + season-correction
  [
    ask patches
    [
      if-else member? (days - season-correction) rainy-days
      [
; In the version for sensitivity analysis, these probabilities are dependent on each other, so that the local minimum is always higher than the slope and the maximum
; differences between topographical objects are the bigger, the less soil, because of its properties, can retain
        if topography = "local maximum" [set p-water-availability (av-water-min - av-water-dif-min-slope - av-water-dif-slope-max) if p-water-availability < 0 [set p-water-availability 0] ]
        if topography = "slope" [set p-water-availability (av-water-min - av-water-dif-min-slope) if p-water-availability < 0 [set p-water-availability 0]]
        if topography = "local minimum" [set p-water-availability av-water-min]
      ]
      [
; If dew matters in the model version, there is some extra water available at the local minima and slopes
        if-else dew? = TRUE
      [
          set p-water-availability (av-water-min - av-water-dif-min-slope - av-water-dif-slope-max - av-water-dif-max-dew)
          if p-water-availability < 0 [set p-water-availability 0]
        ]
        [set p-water-availability 0]
      ]
    ]
; If it is a vegetation season:
  germination    ; Seeds: daily
  competition    ; Plants: daily
  reproduction   ; Plants: under certain conditions
  ]
  [ask patches [set p-water-availability 0]]
end

to dispersal
  ; The FIRST TYPE of dispersal; escape strategy; seeds are dispersed with, for example, wind, randomly, according to species-specific dispersal kernels
  ask plants with [adult? = false and s-disp-type = 1 and s-dispersed? = false and age >= s-disp-age]
    [
      set s-dispersed? true
      rt random-float 360
      ; Weibull-shaped dispersal, parameters β (s-disp-shape) and δ (s-disp-scale) are species-specific
      let s-disp-distance s-disp-scale * ((- ln (1 - random-float 1)) ^ (1 / s-disp-shape))
      jump s-disp-distance * size-factor
    ]
  ; The SECOND TYPE of dispersal, ex. Anastatica; protection strategy; seeds are dispersed on rainy days
ask patches
  [
    if member? days rainy-days
    [
      ask plants-here with [adult? = false and s-disp-type = 2 and s-dispersed? = false and age >= s-disp-age]
      [
        set s-dispersed? true
        rt random-float 360
        let s-disp-distance s-disp-scale * ((- ln (1 - random-float 1)) ^ (1 / s-disp-shape))
        jump s-disp-distance * size-factor
      ]
    ]
  ]
end

to ad-mortality
  ask plants with [adult? = true]
  [
 ; This is a limitation, because of the carrying capacity of a patch, in this case, not more than 10 plants per 1 dm^2
   if count plants-here with [adult? = true] > 1000 [die]
    if ad-biomass < s-mass [die]
  ]
  ask plants with [adult? = true and profile = "Anastatica"]
  [
 ; Some plants die randomly from natural reasons
    if-else ad-biomass >= ad-biomass-max * juv-biomass-an

 ; For size-dependent mortality, plants are divided into two groups: juvenile and adult depending on the fraction (juv-biomass) of the maximum biomass they have as their biomass
 ; if the biomass is less then this fraction (the plant is juvenile), the plant has fewer chances to survive (the probability to survive is corrected by the factor juv-surv)
    [ if random-float 1 >= ad-prob-surv  [die] ]
    [ if random-float 1 >= ad-prob-surv * juv-surv-an [die]]
  ]
  ask plants with [adult? = true and profile = "Malva"]
  [
 ; Some plants die randomly from natural reasons
    if-else ad-biomass >= ad-biomass-max * juv-biomass-ma
 ; For size-dependent mortality, plants are divided into two groups: juvenile and adult depending on the fraction (juv-biomass) of the maximum biomass they have as their biomass
 ; if the biomass is less then this fraction (the plant is juvenile), the plant has fewer chances to survive (the probability to survive is corrected by the factor juv-surv)
    [ if random-float 1 >= ad-prob-surv  [die] ]
    [ if random-float 1 >= ad-prob-surv * juv-surv-ma [die]]
  ]
end


to germination
  ask patches
  [
    ask plants-here with [adult? = false and s-dispersed? = true]
    [
 ; This is a limitation, because of the carrying capacity of a patch, in this case,
 ; not more than 1 seeds per 1cm^2 is allowed
   if count plants-here with [adult? = false] > 10000 [die]
 ; It is checked if there is enough soil water for germination if they've finished their dormancy period and other criteria are met
      if random-float 1  < p-water-availability AND random-float 1 < s-pr-germ AND age >= s-days-dorm
      [
        set adult? true
        set age 1
        if-else zoi-approach
        [
        ; Current biomass is calculated after Weiner & Damgaard, 2006
          set ad-biomass s-mass + ad-rgr * (1 - (s-mass ^ 2 / ad-biomass-max ^ (4 / 3)))
        ]
        [
        ; Current biomass is calculated after Lin et al., 2012
          set ad-biomass s-mass + ad-rgr  * (1 - ((s-mass / ad-biomass-max) ^ (1 / 4)))
        ]
      ]
    ]
  ]
end

to competition
  ask patches
  [
    set p-biomass_acc 0
 ; This line is only for visualization
    set pcolor 39
  ]
 ; Step 1: Plants "occupy" their surrounding
  ask plants with [adult? = true]
  [
    if-else zoi-approach
    [
    ; The radius of ZOI is calculated after Weiner & Damgaard, 2006
      set ad-radius sqrt ( (ad-biomass ^ ( 2 / 3) ) / pi)
    ]
    [
    ; The radius of ZOI is calculated after Lin et al., 2012
      set ad-radius ( ad-biomass ^ ( 3 / 8 ) ) * ( pi ^ ( -1 / 2 ) )
    ]
   ; Calculates the number of patches inside the radius
    set ad-zoi-patches (count patches in-radius (ad-radius * size-factor))
   ; Assigns to the patches in the radius of ZOI biomass accumulation value
    ask patches in-radius (ad-radius * size-factor) [set p-biomass_acc p-biomass_acc + [ad-biomass ^ ad-comp-asymmetry] of myself
    set pcolor pcolor - 1  ; This line is only for visualization
    ]
  ]
; Step 2: Inversion of p-biomass_acc (to speed up the model), after Lin et al., 2012
  ask patches with [p-biomass_acc > 0][set p-biomass_acc_inv 1 / p-biomass_acc]
; Step 3: Plants check how much resources they obtain
  ask plants with [adult? = true]
  [
    if ad-zoi-patches > 0
    [
    ; Share of resources in the radius of ZOI
      let ad-zoi-resource-share (sum ([p-biomass_acc_inv * [ad-biomass ^ ad-comp-asymmetry] of myself] of patches in-radius (ad-radius * size-factor)))
    ; Available resources per occupied patch
      let ad-ratio-resource-area ad-zoi-resource-share / ad-zoi-patches
      if-else zoi-approach
      [
      ; The effective area is calculated after Weiner & Damgaard, 2006
        set ad-zoi-overlap ad-biomass ^ (2 / 3) * ad-ratio-resource-area
      ]
      [
      ; The effective area is calculated after Lin et al., 2012
        set ad-zoi-overlap (ad-biomass ^ (3 / 4) ) * ad-ratio-resource-area
      ]
    ]
  ]
growth
end

to growth
  ask patches
  [
  ; Plants that have the number of patches in their ZOI radius > 0 can grow on the patches
    ask plants-here with [adult? = true and ad-zoi-patches > 0 ]
     [ if p-water-availability > random-float 1
      [
        if-else zoi-approach
        [
        ; Current biomass is calculated after Weiner & Damgaard, 2006
          set ad-biomass ad-biomass + ad-rgr * (ad-zoi-overlap - (ad-biomass ^ 2 / ad-biomass-max ^ (4 / 3)))
        ]
        [
        ; Current biomass is calculated after Lin et al., 2012
          set ad-biomass ad-biomass + ad-rgr * ad-zoi-overlap * (1 - ((ad-biomass / ad-biomass-max) ^ (1 / 4)))
        ]
      ]
    ]
  ]
end

to ageing
  ; Age of both adult plants and seeds is increased with each tick
  ask plants
  [
    set age age + 1
  ]
end

to reproduction
 ; Two conditions are checked: plants are adult and they gained the necessary biomass (can produce at least one seed)
  ask plants with [adult? = true AND count-s-number > 0]
[
    if (age >= ad-age-repr) and (ad-biomass > ad-biomass-repr)
    [
      let parenty ycor
      let parentx xcor
      ; Number of seeds, which can be produced based on current biomass
      set s-number-actual int ( ad-biomass / ad-biomass-repr )
      ; Biomass is reduced depending on the number of seeds produced
      set ad-biomass ad-biomass - s-mass * s-number-actual
      set count-s-number count-s-number - s-number-actual
      if s-number-actual > count-s-number
      [
        set s-number-actual count-s-number
        set count-s-number 0
      ]
      hatch-plants int s-number-actual
      [
        set adult? false
        set ad-biomass s-mass
        set age 1
        set s-dispersed? false
      ]
    ]
    if ad-biomass < s-mass
    ; If as a result of the reproduction process plants lost their biomass below the critical level
    [die]
  ]
end

to count-plants
  set ad-number-An count plants with [profile = "Anastatica" and adult? = true]
  set ad-number-Ma count plants with [profile = "Malva" and adult? = true]
  set s-number-An count plants with [profile = "Anastatica" and adult? = false]
  set s-number-Ma count plants with [profile = "Malva" and adult? = false]

end

to plotting
  ; Calculating daily mean biomass of plants from this specific species
  set-current-plot "Biomass-ma"
  if-else count plants with [profile = "Malva" and adult? = true] >= 1
  [
    set bio-mean-ma mean [ad-biomass] of plants with [profile = "Malva" and adult? = true]
    set-current-plot-pen "biomass-ma" plot bio-mean-ma
  ]
  [set-current-plot-pen "biomass-ma" plot 0]
  set-current-plot "Biomass-an"
  if-else count plants with [profile = "Anastatica" and adult? = true] >= 1
  [
    set bio-mean-an mean [ad-biomass] of plants with [profile = "Anastatica" and adult? = true]
    set-current-plot-pen "biomass-an" plot bio-mean-an
  ]
  [set-current-plot-pen "biomass-an" plot 0]
end
@#$#@#$#@
GRAPHICS-WINDOW
210
10
925
376
-1
-1
7.0
1
10
1
1
1
0
1
1
1
0
100
0
50
0
0
1
ticks
30.0

BUTTON
5
10
60
43
setup
setup-environment
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
66
10
128
43
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
5
100
177
133
in-s-number-an
in-s-number-an
0
10000
9999.0
1
1
NIL
HORIZONTAL

MONITOR
12
292
94
337
s-number-an
s-number-an
17
1
11

PLOT
0
385
200
535
Plants
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"ad-number-an" 1.0 0 -1184463 true "" "plot count plants with [profile = \"Anastatica\" and adult? = true]"
"ad-number-ma" 1.0 0 -14439633 true "" "plot count plants with [profile = \"Malva\" and adult? = true]"

PLOT
240
385
440
535
Seeds
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"s-number-an" 1.0 0 -1184463 true "" "plot count plants with [profile = \"Anastatica\" and adult? = false]"
"s-number-ma" 1.0 0 -14439633 true "" "plot count plants with [profile = \"Malva\" and adult? = false]"

MONITOR
5
46
62
91
NIL
years
17
1
11

SWITCH
930
30
1075
63
dew?
dew?
0
1
-1000

SWITCH
930
63
1075
96
zoi-approach
zoi-approach
1
1
-1000

TEXTBOX
930
100
1090
141
On - Weiner & Damgaard, 2006\nOff - Lin et al., 2012
11
0.0
1

SLIDER
932
231
1077
264
min-height
min-height
0
50
6.0
1
1
NIL
HORIZONTAL

SLIDER
932
264
1077
297
height-dif
height-dif
0
50
9.0
1
1
NIL
HORIZONTAL

SLIDER
932
297
1077
330
sd-height
sd-height
0
1
0.44
0.01
1
NIL
HORIZONTAL

MONITOR
932
330
1077
375
local minimum
count patches with [topography = \"local minimum\"]
17
1
11

MONITOR
933
374
1077
419
local maximum
count patches with [topography = \"local maximum\"]
17
1
11

TEXTBOX
965
214
1041
232
TOPOGRAPHY
11
0.0
1

TEXTBOX
1106
15
1256
33
PRECIPITATION
11
0.0
1

SLIDER
1087
29
1259
62
min-days-year
min-days-year
0
50
9.0
1
1
NIL
HORIZONTAL

SLIDER
1087
62
1259
95
rain-days-year-dif
rain-days-year-dif
0
50
32.0
1
1
NIL
HORIZONTAL

SLIDER
1087
95
1269
128
rain-season
rain-season
1
200
110.0
1
1
NIL
HORIZONTAL

PLOT
460
385
660
535
Biomass-an
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"biomass-an" 1.0 0 -2674135 true "" ""

SLIDER
1100
235
1272
268
av-water-min
av-water-min
0.5
1
0.79
0.01
1
NIL
HORIZONTAL

SLIDER
1100
265
1275
298
av-water-dif-min-slope
av-water-dif-min-slope
0
0.5
0.26
0.01
1
NIL
HORIZONTAL

SLIDER
1100
295
1275
328
av-water-dif-slope-max
av-water-dif-slope-max
0
0.5
0.29
0.01
1
NIL
HORIZONTAL

SLIDER
1100
325
1275
358
av-water-dif-max-dew
av-water-dif-max-dew
0
0.5
0.24
0.01
1
NIL
HORIZONTAL

TEXTBOX
1105
211
1350
246
WATER AVAILABLE FOR PLANTS
11
0.0
1

SLIDER
1090
130
1262
163
dew-days-before
dew-days-before
1
30
9.0
1
1
NIL
HORIZONTAL

SLIDER
1090
165
1262
198
dew-days-after
dew-days-after
1
30
9.0
1
1
NIL
HORIZONTAL

SLIDER
1290
25
1462
58
ad-rgr-an
ad-rgr-an
0.01
1
0.56
0.01
1
NIL
HORIZONTAL

TEXTBOX
1290
10
1440
28
Anastatica traits
11
0.0
1

SLIDER
1290
55
1462
88
ad-biomass-max-an
ad-biomass-max-an
1
1000
466.0
1
1
NIL
HORIZONTAL

SLIDER
1290
85
1462
118
ad-prob-surv-an
ad-prob-surv-an
0
1
0.55
0.05
1
NIL
HORIZONTAL

SLIDER
1290
115
1462
148
s-sp-number-an
s-sp-number-an
1
100
47.0
1
1
NIL
HORIZONTAL

SLIDER
1290
175
1462
208
ad-age-repr-an
ad-age-repr-an
1
150
82.0
1
1
NIL
HORIZONTAL

SLIDER
1290
205
1462
238
s-disp-age-an
s-disp-age-an
0
100
41.0
1
1
NIL
HORIZONTAL

SLIDER
1290
235
1462
268
s-disp-scale-an
s-disp-scale-an
0
10
5.0
1
1
NIL
HORIZONTAL

SLIDER
1290
265
1462
298
s-disp-shape-an
s-disp-shape-an
0.1
1.5
0.55
0.05
1
NIL
HORIZONTAL

SLIDER
1290
295
1462
328
s-days-dorm-an
s-days-dorm-an
1
365
217.0
1
1
NIL
HORIZONTAL

SLIDER
1290
325
1462
358
s-pr-germ-an
s-pr-germ-an
0
1
0.45
0.05
1
NIL
HORIZONTAL

SLIDER
1290
355
1460
388
ad-comp-asymmetry-an
ad-comp-asymmetry-an
0
1
0.5
0.05
1
NIL
HORIZONTAL

SLIDER
1290
385
1462
418
s-disp-type-an
s-disp-type-an
1
2
2.0
1
1
NIL
HORIZONTAL

SLIDER
1290
145
1462
178
s-mass-an
s-mass-an
0.1
10
3.5
0.1
1
NIL
HORIZONTAL

MONITOR
8
246
100
291
ad-number-an
ad-number-an
17
1
11

SLIDER
5
165
177
198
in-ad-number-an
in-ad-number-an
0
100
4.0
1
1
NIL
HORIZONTAL

SLIDER
1290
420
1462
453
juv-biomass-an
juv-biomass-an
0
1
0.67
0.01
1
NIL
HORIZONTAL

SLIDER
1290
455
1462
488
juv-surv-an
juv-surv-an
0
1
0.0
0.01
1
NIL
HORIZONTAL

MONITOR
10
335
102
380
bio-mean-an
bio-mean-an
17
1
11

SLIDER
1495
420
1667
453
juv-biomass-ma
juv-biomass-ma
0
1
0.55
0.01
1
NIL
HORIZONTAL

SLIDER
1495
455
1667
488
juv-surv-ma
juv-surv-ma
0
1
1.0
0.01
1
NIL
HORIZONTAL

SLIDER
5
135
177
168
in-s-number-ma
in-s-number-ma
0
10000
10000.0
1
1
NIL
HORIZONTAL

SLIDER
5
200
177
233
in-ad-number-ma
in-ad-number-ma
0
100
3.0
1
1
NIL
HORIZONTAL

MONITOR
100
245
192
290
ad-number-ma
ad-number-ma
17
1
11

MONITOR
100
290
182
335
s-number-ma
s-number-ma
17
1
11

MONITOR
100
335
197
380
NIL
bio-mean-ma
17
1
11

TEXTBOX
1505
5
1655
23
Malva traits
12
0.0
1

SLIDER
1495
25
1667
58
ad-rgr-ma
ad-rgr-ma
0.01
1
0.39
0.01
1
NIL
HORIZONTAL

SLIDER
1495
55
1692
88
ad-biomass-max-ma
ad-biomass-max-ma
1
1000
424.0
1
1
NIL
HORIZONTAL

SLIDER
1495
85
1667
118
ad-prob-surv-ma
ad-prob-surv-ma
0
1
0.45
0.05
1
NIL
HORIZONTAL

SLIDER
1495
115
1667
148
s-sp-number-ma
s-sp-number-ma
1
100
51.0
1
1
NIL
HORIZONTAL

SLIDER
1495
145
1667
178
s-mass-ma
s-mass-ma
0.1
10
3.7
0.1
1
NIL
HORIZONTAL

SLIDER
1495
175
1667
208
ad-age-repr-ma
ad-age-repr-ma
1
150
80.0
1
1
NIL
HORIZONTAL

SLIDER
1495
205
1667
238
s-disp-age-ma
s-disp-age-ma
0
100
51.0
1
1
NIL
HORIZONTAL

SLIDER
1495
235
1667
268
s-disp-scale-ma
s-disp-scale-ma
0
10
4.0
1
1
NIL
HORIZONTAL

SLIDER
1495
265
1667
298
s-disp-shape-ma
s-disp-shape-ma
0.1
1.5
0.7
0.05
1
NIL
HORIZONTAL

SLIDER
1495
295
1667
328
s-days-dorm-ma
s-days-dorm-ma
1
365
189.0
1
1
NIL
HORIZONTAL

SLIDER
1495
325
1667
358
s-pr-germ-ma
s-pr-germ-ma
0
1
0.55
0.05
1
NIL
HORIZONTAL

SLIDER
1495
355
1670
388
ad-comp-asymmetry-ma
ad-comp-asymmetry-ma
0
1
0.5
0.05
1
NIL
HORIZONTAL

SLIDER
1495
385
1667
418
s-disp-type-ma
s-disp-type-ma
1
2
1.0
1
1
NIL
HORIZONTAL

PLOT
680
385
880
535
Biomass-ma
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"biomass-ma" 1.0 0 -16777216 true "" ""

SLIDER
1290
490
1462
523
ad-biomass-repr-an
ad-biomass-repr-an
1
100
50.0
1
1
NIL
HORIZONTAL

SLIDER
1495
490
1667
523
ad-biomass-repr-ma
ad-biomass-repr-ma
1
100
53.0
1
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

This is the individual-based ATID-model (Anastatica Trait-Based and Individual-based Desert model; atid – future in Hebrew) to investigate, which environmental factors and plant functional traits are important for plant performance and community dynamics in deserts. The model combines trait-based (Violle et al., 2007) and individual-based modelling approaches (Grimm and Railsback, 2005).


## HOW IT WORKS

Agents are adult individual plants and seeds of two generic species that are inspired by Anastatica hierochuntica and Malva parviflora. These two desert annuals represent contrasting dispersal strategies: protective and escape strategies. The model has an in-built opportunity to be extended by an unlimited number of other species with their functional trait combination, depending on the plant community that the model represents. Each species in this model is represented as a set of plant functional traits related to plant-plant interactions and to coping with abiotic stress.
The world is made up of patches representing arid sites with different topography, soil water availability and indirectly salinity with a grain of 1m by 1 m. The spatial extent covers an area of 100 m by 50 m. The tick is one day during the simulated vegetation season. The beginning of each simulated year corresponds to the beginning of the vegetation season. Precipitation input is simulated annually.
During each tick, the two main parts of the model, Vegetation and Environment, are executed. The simulation loop continues until the specified time extent is reached or no adult plant or seed agents exist anymore.

There are two groups of agents: seeds and adult plants.
•	Processes at the level of seeds: 
Seed germination transform seeds into plants if there is enough soil water content. Seed dispersal spreads out seeds. Seed mortality represents the loss of seeds caused by different processes. 
•	Processes at the level of adult plants: 
Adult plants experience natural mortality, inter- and intra-specific competition, growth, ageing and reproduction. Natural mortality is the probability to die because of any natural mortality factors other than direct competition, e.g. disturbance by wild animal activity. Competition occurs via overlapping zones of influence. Growth is represented by an increase in biomass, based on the share of resources that each adult plant gets as a result of the run of the competition procedure. Ageing is applied for both seeds and adult plants. Reproduction leads to the production of seeds under certain conditions.
During each tick that is part of the vegetation season, soil water availability is updated for each patch depending on its topography and precipitation. The model considers topography as the most important environmental factor, and all other environmental factors in the model depend on this one.


## HOW TO USE IT

The sliders in the TOPOGRAPHY group (steepness of the site : min-height, height-dif and  smoothness of the site: sd-height) control the topographical characteristics of the world. The monitors “local minimum” and “local maximum” show how many of each of these topographical objects were simulated by the model.
The sliders in the PRECIPITATION group (minimal and maximal number of rainy days in a year: min-days-year, rain-days-year-dif; length of the vegetation season: rain-season; additional days when dew deposition might affect plant growth and seed germination: dew-days-before, dew-days-after) set parameters for the annual precipitation simulation.
The sliders in WATER AVAILABLE FOR PLANTS control how much water is available for plant growth and seed germination in different topographical objects. These values indirectly simulate soil characteristics and soil ability to retain and provide water for plant growth and seed germination. The larger the difference between water availabilities is for specific topographical objects, the more distinguishable are the spatial preferences of the plants.
Two groups of sliders (Anastatica traits and Malva traits) control plant functional traits of two generic plant species representing two contrasting dispersal strategies named after A. hierochuntica and M. parviflora.
The sliders in-s-number-an and in-s-number-ma control the initial number of seeds and the sliders in-ad-number-an and in-ad-number-ma control the initial number of adult plants of each species.
The switch zoi-approach provides a choice between two different algorithms of calculating zone-of-influence and plant growth. 
The switch dew? controls the influence of precipitation in the form of dew on the plant community.


## THINGS TO NOTICE

Watch desert plant community dynamics in the world display. The darker patches show higher density of vegetation. The higher density vegetation might correspond to environmental conditions that are more favourable for plant growth and seed germination. The monitors Plants and Seeds demonstrate two-species community dynamics as the annual fluctuations in numbers of individuals. The monitors Biomass-an and Biomass-ma show how the mean population biomass of adult plants of the corresponding species changes during the vegetation season. One set of environmental conditions might be favourable for large populations with relatively small plants, other sets of environmental conditions for small populations with large adult plants, etc.


## THINGS TO TRY

Change the environmental factors: topographical and precipitation parameters. Watch how the steepness and the smoothness of the site affects plant community dynamics. Observe the effect of rain season length and the frequency of rain events.  This can provide insights on why this or that species dominate in plant communities across the precipitation gradient or with different topographical characteristics.
Change the water availability for different topographical objects, simulating different soil types. Observe how underlying soil affects community dynamics.
Change the values of plant functional traits. Observe which traits are more favourable for plant growth and survival under given environmental conditions.  Watch the plant community with contrasting versus the same seed dispersal strategies. How do these changes affect community dynamics? Change the juvenile biomass threshold and survival probability. Watch how higher or lower juvenile mortality affects the community.
Run the model with different ZOI scenarios, with plant growth and zone-of-influence calculated after Lin et al. 2012 or Weiner & Damgaard, 2006. Compare how it affects plant growth, competition and survival.
Run the model with and without dew effect. Compare the community dynamics, species abundance and spatial distribution.


## EXTENDING THE MODEL

Input derived from real-world data can replace algorithms simulating precipitation and topography. Adding plant species with species-specific sets of functional traits will extend the model and male it multispecies. Different approaches might be needed to calculate zones-of-influence and plant growth. Wherein, the choice in favour of this or that approach depends on the allometric relations derived from empirical measurements. Similarly, new data may help to elaborate on soil water availability and salinity submodels. To capture intraspecific trait variability, the model should include distributions of plant functional traits.
## NETLOGO FEATURES
The model includes inbuilt topography and precipitation simulators.


## RELATED MODELS

•	Lin, Y., Berger, U., Grimm, V., & Ji, Q. R. (2012). Differences between symmetric and asymmetric facilitation matter: exploring the interplay between modes of positive and negative plant interactions. Journal of Ecology, 100(6), 1482-1491.

•	Radny, J., & Meyer, K. M. (2018). The role of biotic factors during plant establishment in novel communities assessed with an agent-based simulation model. PeerJ, 6, e5342.


## CREDITS AND REFERENCES

The current model is available at: https://github.com/lazakharova/PhDThesis

The model description is submitted as: 
•	Zakharova L., Meyer K.M., Seifan M. (20XX) Combining trait- and individual-based modelling to understand desert plant community dynamics. 

Please cite the NetLogo software as: 
•	Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL. 

Other references:
•	Grimm, V., & Railsback, S. F. (2005). Individual-based modeling and ecology (Vol. 8). Princeton university press.
•	Violle, C., Navas, M. L., Vile, D., Kazakou, E., Fortunel, C., Hummel, I., & Garnier, E. (2007). Let the concept of trait be functional!. Oikos, 116(5), 882-892.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

okr
true
0
Circle -7500403 false true 45 45 210
Circle -7500403 false true 55 55 190

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
1
@#$#@#$#@
