parse_primersearch <- function(file_path) {
  # Split output into chunks for each primer--------------------------------------------------------
  raw_output <- readLines(file_path)
  primer_indexes <- grep("Primer name ", raw_output, fixed=TRUE, value=FALSE)
  primer_chunk_id <- findInterval(seq_along(raw_output), primer_indexes)
  primer_chunks <- vapply(split(raw_output, primer_chunk_id)[-1],
                          paste, character(1), collapse = "\n")
  names(primer_chunks) <- stringr::str_match(primer_chunks, "Primer name ([^\n]*)")[,2]
  # Extract amplicon data from each chunk and combine ----------------------------------------------
  pattern <- paste("Amplimer ([0-9]+)",
                   "\tSequence: ([^\n]*)",
                   "\t([^\n]*)",
                   "\t([^\n]+) hits forward strand at ([0-9]+) with ([0-9]+) mismatches",
                   "\t([^\n]+) hits reverse strand at \\[([0-9]+)\\] with ([0-9]+) mismatches",
                   "\tAmplimer length: ([0-9]+) bp", sep = '\n')
  primer_data <- stringr::str_match_all(primer_chunks, pattern)
  if (length(rep(names(primer_chunks), vapply(primer_data, nrow, numeric(1)))) ==1) {
    primer_data <- as.data.frame(cbind(t(rep(names(primer_chunks), vapply(primer_data, nrow, numeric(1)))),
                                       t(do.call(rbind, primer_data)[, -1])), stringsAsFactors = FALSE)
    # Reformat amplicon data -------------------------------------------------------------------------
    colnames(primer_data) <- c("pair_name", "amplimer", "seq_id", "name", "f_primer", "f_index",
                               "f_mismatch",  "r_primer", "r_index", "r_mismatch", "length")
    primer_data <- primer_data[, c("name", "seq_id", "pair_name", "amplimer", "length", 
                                   "f_primer", "f_index", "f_mismatch",
                                   "r_primer", "r_index", "r_mismatch")]
    numeric_cols <- c("amplimer", "length","f_index", "f_mismatch",
                      "r_index", "r_mismatch")
    
    for (col in numeric_cols) primer_data[, col] <- as.numeric(primer_data[, col])  
    return(primer_data)
  }  else {
    primer_data <- as.data.frame(cbind(rep(names(primer_chunks), vapply(primer_data, nrow, numeric(1))),
                                       do.call(rbind, primer_data)[, -1]), stringsAsFactors = FALSE)
    # Reformat amplicon data -------------------------------------------------------------------------
    colnames(primer_data) <- c("pair_name", "amplimer", "seq_id", "name", "f_primer", "f_index",
                               "f_mismatch",  "r_primer", "r_index", "r_mismatch", "length")
    primer_data <- primer_data[, c("name", "seq_id", "pair_name", "amplimer", "length", 
                                   "f_primer", "f_index", "f_mismatch",
                                   "r_primer", "r_index", "r_mismatch")]
    numeric_cols <- c("amplimer", "length","f_index", "f_mismatch",
                      "r_index", "r_mismatch")
    
    for (col in numeric_cols) primer_data[, col] <- as.numeric(primer_data[, col]) 
    return(primer_data)
  }
}   
