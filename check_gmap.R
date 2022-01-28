# some functions from STITCH that I used to debug the genetic map.

fill_in_genetic_map_rate_column <- function(genetic_map) {
    genetic_map[, "COMBINED_rate.cM.Mb."] <- 0
    for(iRow in 1:(nrow(genetic_map) - 1)) {
        distance <-
            genetic_map[iRow + 1, "position"] -             
            genetic_map[iRow, "position"]
        difference <-
            genetic_map[iRow + 1, "Genetic_Map.cM."] -
            genetic_map[iRow, "Genetic_Map.cM."]
        genetic_map[iRow, "COMBINED_rate.cM.Mb."] <-
            difference / (distance / 1e6)
    }
    return(genetic_map)
}

fill_in_genetic_map_cm_column <- function(genetic_map) {
    genetic_map[1, "Genetic_Map.cM."] <- 0
    for(iRow in 2:nrow(genetic_map)) {
        distance <-
            genetic_map[iRow, "position"] -
            genetic_map[iRow - 1, "position"]
        rate <- genetic_map[iRow - 1, "COMBINED_rate.cM.Mb."]
        genetic_map[iRow, "Genetic_Map.cM."] <-
            genetic_map[iRow - 1, "Genetic_Map.cM."] + 
            (distance / 1e6) * rate
    }
    return(genetic_map)
}

validate_genetic_map <- function(genetic_map, verbose = TRUE) {
    if (ncol(genetic_map) != 3) {
        stop(paste0("Provided reference genetic map has ", ncol(genetic_map), " columns, where 3 are expected, with columns that provide position (in bp), combined genetic rate map (in cM/Mbp), and genetic map (in cM)"))
    }
    original_colnames_genetic_map <- colnames(genetic_map)
    colnames(genetic_map) <- c("position", "COMBINED_rate.cM.Mb.", "Genetic_Map.cM.")
    ## do not allow any NAs or non-numeric characters
    if (sum(is.na(genetic_map) | is.na(is.na(genetic_map * 2))) > 0) {
        stop(paste0("Provided reference genetic map contains non-numeric entries or NA entries. Please ensure all entries are numeric and non-empty"))
    }
    ## 
    check_col_3 <- fill_in_genetic_map_cm_column(genetic_map)[, 3]
    x <- (genetic_map[, 3] - check_col_3) > 1e-6
    if (sum(x) > 0) {
        if (verbose) {
            print(genetic_map[c(max(1, x - 1):min(x + 1, nrow(genetic_map))) , ])
        }
        stop(paste0("Error interpreting genetic map around row ", which.max(x), ". Provided genetic map does not satisfy expectations. See above print. Let g be the genetic map with 3 columns. In 1-based coordinates, recall that g[iRow, 3] must equal g[iRow - 1, 3] + ( (g[iRow, 1] - g[iRow - 1, 1]) * g[iRow - 1, 2])"))
    }
    return(NULL)
}

