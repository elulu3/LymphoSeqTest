#' Export sequences in fasta format
#' 
#' Export junction or amino acid sequences in fasta format.
#' 
#' @param sample_table A tibble consisting of antigen receptor sequences 
#' imported by the LymphoSeq function readImmunoSeq.
#' @param type A character vector indicating whether "junction_aa" or "junction" sequences
#' should be exported.  If "junction_aa" is specified, then run productiveSeqs first.
#' @param names A character vector of one or more column names to name the sequences.
#' If "rank" is specified, then the rank order of the sequences by frequency is used.
#' @return Exports fasta files to the working directory.
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' 
#' stable <- readImmunoSeq(path = file_path)
#' 
#' exportFasta(study_table = stable, type = "junction", names = c("junction_aa", "duplicate_count"))
#' 
#' atable <- productiveSeq(study_table = stable, aggregate = "junction_aa")
#' 
#' exportFasta(study_table = atable, type = "junction_aa", names = "duplicate_frequency")
#' @export
#' @import tidyverse
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings AAStringSet
#' @importFrom Biostrings writeXStringSet
exportFasta <- function(study_table, type = "junction", 
                        names = c("rank", "junction_aa", "duplicate_count")) {
    if (type == "junction") {
        study_table <- study_table %>% 
                       arrange(repertoire_id, desc(duplicate_frequency)) %>% 
                       tibble::rowid_to_column() %>% 
                    #    rename(rowid = rank) %>%
                       mutate(sequences = junction) %>%
                       tidyr::unite(fasta_name, names)
    } else if (type == "junction_aa") {
        study_table <- productiveSeq(study_table)
        study_table <- study_table %>% 
                       arrange(repertoire_id, desc(duplicate_frequency)) %>% 
                       tibble::rowid_to_column() %>% 
                    #    rename(rowid = rank) %>%
                       mutate(sequences = junction_aa) %>%
                       tidyr::unite(fasta_name, names) 
    }
    study_table %>% 
    group_by(repertoire_id) %>% 
    group_split() %>% 
    purrr::map(~writeFasta(.x, type))
    message(paste("Fasta files saved to", getwd()))
}

writeFasta <- function(sample_table, type) {
    repertoire_id <- sample_table$repertoire_id[1]
    if (type == "junction"){
        fasta <- Biostrings::DNAStringSet(sample_table$sequences)
    } else if (type == "junction_aa") {
        fasta <- Biostrings::AAStringSet(sample_table$sequences)
    }
    names(fasta) <- sample_table$fasta_name
    Biostrings::writeXStringSet(fasta, paste(repertoire_id, "fasta", sep="."))
}
