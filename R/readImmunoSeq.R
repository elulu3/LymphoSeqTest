#' Read ImmunoSeq files
#'
#' Imports tab-separated value (.tsv) files exported by the Adaptive 
#' Biotechnologies ImmunoSEQ analyzer, BGI IR-SEQ, MiXCR and stores 
#' them as MiAIRR compliant tibble. 
#' 
#' @md
#' @name readImmunoSeq
#' @describeIn readImmunoSeq Read a set of files
#' @param path Path to the directory containing tab-delimited files.  Only
#' files with the extension .tsv are imported.  The names of the data frames are 
#' the same as names of the files.
#' 
#' @return Returns a tibble with MiAIRR headers and repertoire_id
#'
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' 
#' study_table <- readImmunoSeq(path = file.path,
#'                              recursive = FALSE)
#' @export
#' @import tidyverse dplyr stringr
#' @rdname readImmunoSeq
readImmunoSeq <- function(path, recursive = FALSE) {
    if (length(path) > 1) {
        file_paths <- path
    } else if (file_test("-d", path)) {
        file_paths <- list.files(path, 
                                 full.names = TRUE, 
                                 all.files = FALSE, 
                                 recursive = recursive, 
                                 pattern = ".tsv|.txt|.tsv.gz", 
                                 include.dirs = FALSE)
    } else {
        file_paths <- c(path)
    }
    file_num <- length(file_paths)
    file_info <- file.info(file_paths)
    file_paths <- rownames(file_info)[file_info$size > 0]
    if(file_num != length(file_paths)){
        warning("One or more of the files you are trying to import has no sequences and will be ignored.", 
                call. = FALSE)
    }
    airr_headers_path <- system.file("extdata", "AIRR_fields.csv", package = "LymphoSeq2")
    airr_fields <- readr::read_csv(airr_headers_path, trim_ws = TRUE)
    matching_fields <- c(amino_acid = "sequence_aa", aminoAcid = "sequence_aa",
                        aminoAcid.CDR3.in.lowercase. = "sequence_aa", cdr1_amino_acid = "cdr1_aa",
                        cdr1_rearrangement = "cdr1", cdr2_amino_acid = "cdr2_aa",
                        cdr2_rearrangement = "cdr2", cdr3_amino_acid = "cdr3_aa",
                        cdr3_rearrangement = "cdr3", d_allele_ties = "d2_call",
                        dGeneNameTies = "d2_call", count = "duplicate_count",
                        "count (reads)" = "duplicate_count", "count (templates)" = "duplicate_count",
                        "count (templates/reads)" = "duplicate_count", d_resolved = "d_call",
                        dMaxResolved = "d_call", frame_type = "productive", fuction = "productive",
                        j_resolved = "j_call", jMaxResolved = "j_call", locus = "locus",
                        nucleotide = "sequence", nucleotide.CDR3.in.lowercase. = "sequence",
                        v_resolved = "v_call", vMaxResolved = "v_call")
    progress_bar <- progress::progress_bar$new(format = "Reading AIRR-Seq files [:bar] :percent eta: :eta",
                                           total = length(file_paths), clear = FALSE, width = 60)
    progress_bar$tick(0)
    file_list <- file_paths %>% 
                 purrr::map(~ readFiles(.x, airr_fields, matching_fields, progress_bar)) %>%
                 dplyr::bind_rows()
    progress_bar$terminate()
    return(file_list)
}

#' @section Get standard headers:
#'
#' Retrives MiAIRR standard compliant fields from the clone files
#'
#' @param clone_file A .tsv file to read in and standardize its fields to be MiAIRR compliant
#' @param airr_fields A character vector of MiAIRR headers
#' @return Tibble of given data with MiAIRR fields
#'
#' @import tidyverse
#' @export
#' @rdname readImmunoSeq
getStandard <- function(clone_file, airr_fields, matching_fields) {
    clone_data <- fread(clone_file, sep = '\t', select = c(airr_fields, names(matching_fields)),
                        na.strings = c("", "NA", "Nan", "NaN", "unresolved"),
                        showProgress = FALSE, verbose = FALSE) %>%
                    tibble::as_tibble()
    existing_match <- airr_fields[airr_fields %in% colnames(clone_data)]
    if (length(existing_match) == 144) {
        return(clone_data)
    }
    existing_airr_data <- clone_data %>%
                            dplyr::select(existing_match)
    match <- matching_fields[colnames(clone_data)] %>%
             stats::na.omit()
    renamed_data <- clone_data[names(match)] %>%
                dplyr::rename(!!! stats::setNames(names(match), match))
    if ("d2_call" %in% colnames(renamed_data)) {
        renamed_data <- renamed_data %>% 
                        tidyr::separate("d2_call", c("d_call_2", "d2_call"), ",") %>% 
                        dplyr::mutate(d_call = coalesce(d_call, d_call_2)) %>%
                        dplyr::select(-d_call_2)
    }
    clone_data <- dplyr::bind_cols(existing_airr_data, renamed_data)

    file_name <- basename(clone_file)
    file_name <- file_name %>% stringr::str_sub(1, stringr::str_length(file_name) - 4)
    
    clone_data <- clone_data %>% mutate(repertoire_id = file_name)
    if ("sequence" %in% colnames(clone_data) && "sequence_aa" %in% colnames(clone_data)) {
        clone_data <- clone_data %>%
                    dplyr::mutate(junction = sequence,
                                  junction_aa = sequence_aa)
    } else {
        clone_data <- clone_data %>%
                    dplyr::mutate(sequence = junction,
                                  sequence_aa = junction_aa)
    }
    clone_data <- clone_data %>%
                    dplyr::mutate(clone_id = junction)
    unique_nucl <- clone_data %>%
                    dplyr::select(clone_id) %>%
                    dplyr::distinct()
    nucl <- unique_nucl[, "clone_id", drop = TRUE]
    clone_data <- clone_data %>%
            dplyr::mutate(clone_id = stringr::str_c(file_name, "_", c(which(nucl == clone_id)))) 

    return(clone_data)
}

#' @section Read clone file from path:
#'
#' Given the path to a single AIRRSeq clone file, generate MiAIRR compliant tibble.
#'
#' @param clone_file .tsv file containing results from AIRRSeq pipeline
#'
#' @return Tibble in MiAIRR format
#'
#' @import tidyverse
#' @export
#' @rdname readImmunoSeq
readFiles <- function(clone_file, empty_airr_frame, matching_fields, progress) {
    progress$tick()
    file_info <- getStandard(clone_file, colnames(empty_airr_frame), matching_fields)
    if (ncol(file_info) == 144) {
        return(file_info)
    }
    clone_frame <- dplyr::right_join(empty_airr_frame, file_info)
    options(readr.show_progress = FALSE)
        clone_frame <- clone_frame %>%
                    dplyr::mutate(sequence_id = row_number(),
                        junction_length = stringr::str_length(junction),
                        junction_aa_length = stringr::str_length(junction_aa),
                        rev_comp = FALSE,
                        stop_codon = dplyr::if_else(stringr::str_detect(sequence, "\\*") | stringr::str_detect(sequence_aa, "\\*") | 
                                                is.na(sequence) | is.na(sequence_aa), TRUE, FALSE),
                        productive = dplyr::if_else(stringr::str_detect(sequence, "\\*") | stringr::str_detect(sequence_aa, "\\*") | 
                                                is.na(sequence) | is.na(sequence_aa), FALSE, TRUE),
                        v_call = dplyr::if_else(stringr::str_detect(v_call, "/"), "unresolved", v_call),
                        j_call = dplyr::if_else(stringr::str_detect(j_call, "/"), "unresolved", j_call),
                        complete_vdj = dplyr::if_else(is.na(v_call) | is.na(d_call) | is.na(j_call) | 
                                                stringr::str_detect(v_call, "unresolved") | stringr::str_detect(d_call, "unresolved") |
                                                stringr::str_detect(j_call, "unresolved"), FALSE, TRUE),
                        duplicate_frequency = duplicate_count / sum(duplicate_count),
                        reading_frame = dplyr::if_else(stringr::str_detect(sequence_aa, "\\*"), "out-of-frame", "in-frame"),
                        v_family = dplyr::if_else(stringr::str_detect(v_call, "(TRB|TCRB)V"),
                                                stringr::str_extract(v_call, "(TRB|TCRB)V\\d+"), "unrecognized"),
                        j_family = dplyr::if_else(stringr::str_detect(j_call, "(TRB|TCRB)J"),
                                                stringr::str_extract(j_call, "(TRB|TCRB)J\\d+"), "unrecognized"),
                        d_family = dplyr::if_else(stringr::str_detect(d_call, "(TRB|TCRB)D"),
                                                stringr::str_extract(d_call, "(TRB|TCRB)D\\d+"), "unrecognized"))

    return(clone_frame)
}


#' @section Get iReceptor standard format:
#'
#' Returns a tibble that is compliant with the iReceptor repertoire format.
#'
#' @param clone_frame An AIRR compliant tibble
#' @return Tibble in iReceptor format
#'
#' @import tidyverse
#' @export
#' @rdname readImmunoSeq
iReceptorFormat <- function(clone_frame) {
    iReceptor_frame <- clone_frame %>%
                        dplyr::select(-c(repertoire_id, sample_processing_id,
                            data_processing_id, reading_frame, v_family,
                            d_family, j_family, duplicate_frequency))
    return(iReceptor_frame)
}
