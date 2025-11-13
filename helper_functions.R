# Helper functions ####
# for extracting site names from package names
extract_site_names <- function(text){
  str_remove_all(text, "nhmp|hedgehog|[:digit:]|[:punct:]|[:space:]")
}
# for extracting subsite names from deployments$locationNames
extract_subsite_names <- function(text){
  str_remove_all(text, "[:digit:]|[:space:]") %>%
    str_split("_") %>%
    map_chr(\(x) tail(x, 1)) %>%
    str_remove_all("[:punct:]")
}
# for splitting datapackages
split_camtrapDP <- function(package, by){
  gps <- package$data$deployments %>% 
    dplyr::reframe(unique({{by}})) %>%
    dplyr::pull()
  res <- gps %>%
    purrr::map(\(gp) subset_deployments(package, {{by}}==gp))
  names(res) <- gps
  res
}

# add/modify fields in deployments table, 
# and optionally add them to observations table
mutate_camtrapDP <- function(package, ..., obs=TRUE){
  package$data$deployments <- package$data$deployments %>%
    dplyr::mutate(...)
  if(obs){
    newFields <- names(quos(...))
    deps <- dplyr::select(package$data$deployments, 
                          all_of(c("deploymentID", newFields)))
    package$data$observations <- package$data$observations %>%
      dplyr::left_join(deps, by="deploymentID")
  }
  package
}

merge_camtrapDP <- function(pkgs){
  profile <- map(pkgs, \(pkg) pkg$profile)
  name <- map(pkgs, \(pkg) pkg$name)
  id <- map(pkgs, \(pkg) pkg$id)
  created <- map(pkgs, \(pkg) pkg$created)
  image <- map(pkgs, \(pkg) pkg$image)
  sources <- map(pkgs, \(pkg) pkg$sources)
  project <- map(pkgs, \(pkg) pkg$project)
  spatial <- map(pkgs, \(pkg) pkg$spatial)
  directory <- map(pkgs, \(pkg) pkg$directory)
  
  contributors <- map(pkgs, \(pkg) pkg$contributors) %>%
    unlist(recursive = FALSE)
  ids <- map_chr(contributors, \(cn) paste(cn, collapse=" "))
  contributors <- contributors[!duplicated(ids)]
  names(contributors) <- NULL
  
  temporal_start <- map(pkgs, \(pkg) pkg$temporal$start) %>%
    unlist() %>%
    lubridate::as_date() %>%
    min()
  temporal_end <- map(pkgs, \(pkg) pkg$temporal$end) %>%
    unlist() %>%
    lubridate::as_date() %>%
    max()
  temporal <- list(start = temporal_start,
                   end = temporal_end)
  
  taxonomic <- map(pkgs, \(pkg) pkg$taxonomic) %>%
    unlist(recursive = FALSE)
  ids <- map_chr(taxonomic, \(tx) tx$taxonID)
  taxonomic <- taxonomic[!duplicated(ids)]
  names(taxonomic) <- NULL
  
  resources <- map(pkgs, \(pkg) pkg$resources) %>%
    unlist(recursive = FALSE)
  nms <- map_chr(resources, \(re) re$name)
  resources <- resources[!duplicated(nms)]
  names(resources) <- NULL
  
  deps <- lapply(pkgs, function(pkg) pkg$data$deployments)
  meds <- lapply(pkgs, function(pkg) pkg$data$media)
  obss <- lapply(pkgs, function(pkg) pkg$data$observations)
  poss <- lapply(pkgs, function(pkg) pkg$data$positions)
  
  nrow2 <- function(x) if(is.null(x)) 0 else nrow(x)
  ndeps <- unlist(lapply(deps, nrow))
  nmeds <- unlist(lapply(meds, nrow2))
  nobss <- unlist(lapply(obss, nrow))
  nposs <- unlist(lapply(poss, nrow2))
  
  f <- function(d) if(nrow(d)==0) NULL else d
  meds <- lapply(meds, f)
  obss <- lapply(obss, f)
  poss <- lapply(poss, f)
  
  ids <- unlist(lapply(pkgs, function(pkg) pkg$id))
  names <- unlist(lapply(pkgs, function(pkg) pkg$name))
  
  deps <- deps %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(packageName = rep(names, ndeps),
                  packageID = rep(ids, ndeps))
  meds <- meds %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(packageName = rep(names, nmeds),
                  packageID = rep(ids, nmeds))
  obss <- obss %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(packageName = rep(names, nobss),
                  packageID = rep(ids, nobss))
  poss <- poss %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(packageName = rep(names, nposs),
                  packageID = rep(ids, nposs))
  
  data <- list(deployments = deps,
               media = meds,
               observations = obss,
               positions = poss)
  names(pkgs$Dorset 1)
  
  list(profile=profile, name=name, id=id, created=created, 
       contributors=contributors, image=image, cources=sources, 
       project=project, spatial=spatial, temporal=temporal, 
       taxonomic=taxonomic, resources=resources, directory=directory, 
       data=data)
}

add_species <- function(to, from, species){
  toDeps <- to$data$deployments
  depLookup <- from$data$deployments %>%
    select(deploymentID, locationName) %>%
    mutate(deploymentID2 = toDeps$deploymentID[match(locationName,
                                                     toDeps$locationName)])
  obs <- from$data$observations %>%
    dplyr::filter(scientificName %in% species) %>%
    mutate(deploymentID = depLookup$deploymentID2[match(deploymentID,
                                                        depLookup$deploymentID)])
  pos <- from$data$positions %>%
    dplyr::filter(scientificName %in% species) %>%
    mutate(deploymentID = depLookup$deploymentID2[match(deploymentID,
                                                        depLookup$deploymentID)])
  
  to$data$observations <- to$data$observations %>%
    bind_rows(obs)
  to$data$positions <- to$data$positions %>%
    bind_rows(pos)
  
  sppTo <- map_chr(to$taxonomic, \(x) x$scientificName)
  sppFrom <- map_chr(from$taxonomic, \(x) x$scientificName)
  spp <- species[!species %in% sppTo]
  to$taxonomic <- c(to$taxonomic, from$taxonomic[match(spp, sppFrom)])

  to
}

subset_observations <- function(package, choice, suffix = ""){
  package$name <- paste(package$name, suffix, sep = "-")
  package$data$observations <- dplyr::filter(package$data$observations, {{choice}})
  package$data$positions <- dplyr::filter(package$data$positions, {{choice}})
  package
}

replace_year <- function(package) {
  # Leaving for reference, the deployments are updated from the effort
  # spreadsheet
  # MR: I reinstated following 2 lines
  package$data$observations$timestamp[year(package$data$observations$timestamp) < 2023] <- 
    update(package$data$observations$timestamp[year(package$data$observations$timestamp) < 2023], year = 2023)
  
  package$data$deployments$start[year(package$data$deployments$start) < 2023] <- 
    update(package$data$deployments$start[year(package$data$deployments$start) < 2023], year = 2023)
  
  package$data$deployments$end[year(package$data$deployments$end) < 2023] <- 
    update(package$data$deployments$end[year(package$data$deployments$end) < 2023], year = 2023)
  
  package
}

replace_month <- function(package, from, to) {
  # Leaving for reference, the deployments are updated from the effort
  # spreadsheet
  # MR: I reinstated following 2 lines
  package$data$observations$timestamp[month(package$data$observations$timestamp) == from] <- 
    update(package$data$observations$timestamp[month(package$data$observations$timestamp) == from], month = to)
  
  package$data$deployments$start[month(package$data$deployments$start) == from] <- 
    update(package$data$deployments$start[month(package$data$deployments$start) == from], month = to)
  
  package$data$deployments$end[month(package$data$deployments$end) == from] <- 
    update(package$data$deployments$end[month(package$data$deployments$end) == from], month = to)
  
  package
}

switch_am_pm <- function(pkg, earliest="06:00:00", latest="20:00:00"){
  earliest <- lubridate::hms(earliest)
  latest <- lubridate::hms(latest)
  times <- lubridate::hms(format(pkg$data$deployments$start, "%H:%M:%S"))
  depIDs <- pkg$data$deployments$deploymentID
  diffs <- dplyr::case_when(times < earliest ~ 12*60^2,
                            times > latest ~ -12*60^2,
                            .default = 0)
  # Leaving for reference, the deployments are updated from the effort
  # spreadsheet
  #
  #pkg$data$deployments <- pkg$data$deployments %>%
  #  mutate(start = start + diffs,
  #         end = end + diffs)
  pkg$data$media <- pkg$data$media %>%
    mutate(timestamp = timestamp + diffs[match(deploymentID, depIDs)])
  pkg$data$observations <- pkg$data$observations %>%
    mutate(timestamp = timestamp + diffs[match(deploymentID, depIDs)])
  pkg
}

# camtrapDensity correct_time
correct_time2 <- function(package, depID=NULL, locName=NULL, wrongTime, rightTime){
  nullsum <- sum(is.null(depID), is.null(locName))
  if(nullsum!=1)
    stop("One but not both of depID and locName must be provided")
  
  if(!is.null(depID)){
    if(length(depID) != length(wrongTime) | length(depID) != length(rightTime))
      stop("depID, wrongTime and rightTime must all have the same length")
    if(!all(depID %in% package$data$deployments$deploymentID))
      stop("Can't find all depID in package$data$deployments$deploymentID")
  }
  
  if(!is.null(locName)){
    if(length(locName) != length(wrongTime) | length(locName) != length(rightTime))
      stop("locName, wrongTime and rightTime must all have the same length")
    if(!all(locName %in% package$data$deployments$locationName))
      stop("Can't find all locName in package$data$deployments$locationName")
    depID <- with(package$data$deployments, 
                  deploymentID[match(locName, locationName)])
    if(length(depID) != length(locName))
      stop(paste("One or more locations have multiple deployments, use depID instead of locName"))
  }
  
  td <- c(0, difftime(rightTime, wrongTime, tz="UTC"))
  iDep <- 1 + match(package$data$deployments$deploymentID, depID, nomatch = 0)
  package$data$deployments <- package$data$deployments %>%
    dplyr::mutate(start = start + td[iDep],
                  end = end + td[iDep])
  
  iObs <- 1 + match(package$data$observations$deploymentID, depID, nomatch = 0)
  package$data$observations <- package$data$observations %>%
    dplyr::mutate(timestamp = timestamp + td[iObs])
                  
  if("media" %in% names(package$data)){
    iMed <- 1 + match(package$data$media$deploymentID, depID, nomatch = 0)
    package$data$media <- package$data$media %>%
      dplyr::mutate(timestamp = timestamp + td[iMed])
  }
  
  package$temporal$start <- lubridate::as_date(min(package$data$deployments$start))
  package$temporal$end <- lubridate::as_date(max(package$data$deployments$end))
  package
}

# add species from one package to another
# assuming from the same survey but processed separately
# NB media are not added so mediaID values not matched for added observations
add_species <- function(to, from, species){
  fromSpp <- unlist(map(from$taxonomic, \(sp) sp$scientificName))
  toSpp <- unlist(map(to$taxonomic, \(sp) sp$scientificName))
  i <- match(species[!species %in% toSpp], fromSpp)
  if(length(i) > 0)
    to$taxonomic <- c(to$taxonomic, from$taxonomic[i])
  
  toDeps <- to$data$deployments
  depLookup <- from$data$deployments %>%
    select(deploymentID, locationName) %>%
    mutate(deploymentID2 = toDeps$deploymentID[match(locationName,
                                                     toDeps$locationName)])
  obs <- from$data$observations %>%
    dplyr::filter(scientificName %in% species) %>%
    mutate(deploymentID = depLookup$deploymentID2[match(deploymentID,
                                                        depLookup$deploymentID)])
  pos <- from$data$positions %>%
    dplyr::filter(scientificName %in% species) %>%
    mutate(deploymentID = depLookup$deploymentID2[match(deploymentID,
                                                        depLookup$deploymentID)])
  
  to$data$observations <- to$data$observations %>%
    bind_rows(obs)
  to$data$positions <- to$data$positions %>%
    bind_rows(pos)
  
  to
}
