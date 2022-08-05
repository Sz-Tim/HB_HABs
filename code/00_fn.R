## -----------------------------------------------------------------------------
## HAB forecasts
## Functions
## Tim Szewczyk
## -----------------------------------------------------------------------------


get_os <- function() {
  # https://www.r-bloggers.com/2015/06/identifying-the-os-from-r/
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}


setPartTrackDirs <- function(base.dir="./", 
                             java.dir="/home/sa04ts/particle_tracking/") {
  
}

setPartTrackProperties <- function(
  destinationDirectory="out/",
  datadir="/media/archiver/common/sa01da-work/WeStCOMS2/Archive/",
  datadirPrefix="netcdf_",
  datadirSuffix="",
  datadir2="",
  datadir2Prefix="",
  datadir2Suffix="",
  mesh1="/home/sa04ts/FVCOM_meshes/WeStCOMS2_mesh.nc",
  mesh1Type="FVCOM",
  mesh2="",
  coordRef="OSGB1936",
  location="minch",
  minchVersion=2,
  habitat="",
  suffix="",
  sitefile="large_elem_centroids_v2.dat",
  sitefileEnd="fsa_sites_v2.dat",
  verboseSetUp="true",
  start_ymd=20190401,
  end_ymd=20220601,
  numberOfDays=0,
  backwards="false",
  checkOpenBoundaries="false",
  readHydroVelocityOnly="false",
  duplicateLastDay="true",
  recordsPerFile1=25,
  dt=3600,
  verticalDynamics="true",
  fixDepth="false",
  maxDepth="",
  parallelThreads=4,
  restartParticles="",
  restartParticlesCutoffDays=21,
  releaseScenario=1,
  nparts=1,
  setStartDepth="true",
  startDepth=1,
  seasonalDensityPath="",
  thresh=500,
  endOnArrival="false",
  rk4="true",
  stepsPerStep=30,
  diffusion="true",
  variableDiffusion="true",
  D_h=0.1,
  D_hVert=0.001,
  salinityThreshold=0,
  mortalityRate=0,
  salinityMort="false",
  swimLightLevel="false",
  vertSwimSpeedMean=0,
  vertSwimSpeedStd=0,
  sinkingRateMean=0,
  sinkingRateStd=0,
  viabletime=-1,
  maxParticleAge=168,
  viableDegreeDays=-1,
  maxDegreeDays=-1,
  recordPsteps="false",
  splitPsteps="true",
  pstepsInterval=168,
  recordMovement="true",
  recordElemActivity="false",
  recordConnectivity="true",
  connectivityInterval=24,
  recordLocations="true",
  recordArrivals="true"
) {
  args <- formals()
  return(paste(names(args), args, sep="=", collapse="\n"))
}




clean_fsa <- function(x, v1_end) {
  x %>%
    as_tibble %>% 
    filter(!is.na(date_collected)) %>%
    filter(karenia_mikimotoi >= 0) %>% # -99 in karenia...?
    mutate(datetime_collected=as_datetime(date_collected)) %>%
    filter(datetime_collected >= "2013-07-20") %>% 
    mutate(date=date(datetime_collected),
           hour=pmin(20, pmax(5, hour(datetime_collected))),
           grid=if_else(date_collected < v1_end, 1, 2),
           site.id=as.numeric(factor(paste(sin, area, site)))) %>%
    rename(obs.id=oid) %>%
    group_by(sin, area, site) %>% 
    mutate(lon=median(easting), lat=median(northing)) %>%
    ungroup %>%
    mutate(depth_recorded=depth %>% 
             str_to_lower() %>% 
             str_remove_all("m| |\\r|\\n|j|<|>|\\+") %>%
             str_replace("su?.face", "0") %>%
             str_remove("\\(0\\)") %>%
             str_replace("bucket", "NA") %>%
             str_replace("30c", ".3") %>%
             str_replace("55", "5.5") %>%
             str_remove("[0-9]-") %>%
             str_replace("-", "NA") %>%
             as.numeric) %>%
    select(-geom, -easting, -northing, -tide, -datetime_collected, -depth)
}