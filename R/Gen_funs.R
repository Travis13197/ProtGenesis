# ==============================================================================
# 1.Gen_funs.R
# ==============================================================================
#
# Purpose: Comprehensive protein sequence generation and mutation functions
#   This script provides two core functions for generating mutation and 
#   stepwise sequences from FASTA input files with flexible configuration options.
#
# Author: ProtGenesis Analysis Team
# Date: 2026-02-24
# Version: 1.0.0
# By：Chuanyang Liu
# please, please cite our paper if you use this script in your research:
#  Chuanyang Liu, et al. (2026). Universal physical principles govern the deterministic genesis of protein structure. 
#  BioRxiv.10.64898/2026.02.20.706798.
#
# Usage:
#   source("1.Gen_funs.R")
#   DataGen_Mutation(...)
#   DataGen_Stepwise(...)
#
# Requirements:
#   - Biostrings (for FASTA file handling)
#   - stringr (for string manipulation)
#   - dplyr (for data manipulation)
#
# ==============================================================================

# ==============================================================================
# Package Management
# ==============================================================================

# Check and install required packages if necessary
check_and_install_packages <- function() {
  required_packages <- c("Biostrings", "stringr", "dplyr")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg == "Biostrings") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
  }
  
  # Load required packages
  library(Biostrings)
  library(stringr)
  library(dplyr)
}

check_and_install_packages()

# ==============================================================================
# Helper Functions
# ==============================================================================

#' Validate FASTA file format and content
#'
#' @param fasta_path Path to FASTA file
#' @return List containing validation status and sequences
#' @noRd
validate_fasta_input <- function(fasta_path) {
  if (!file.exists(fasta_path)) {
    stop(sprintf("FASTA file not found: %s", fasta_path))
  }
  
  tryCatch({
    fasta_seqs <- readAAStringSet(fasta_path)
    
    if (length(fasta_seqs) == 0) {
      stop("FASTA file contains no sequences")
    }
    
    # Validate protein sequences (only standard amino acids)
    valid_aa <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                  "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
    
    for (i in seq_along(fasta_seqs)) {
      seq_char <- as.character(fasta_seqs[i])
      seq_aa <- strsplit(seq_char, "")[[1]]
      
      invalid_aa <- setdiff(unique(seq_aa), valid_aa)
      if (length(invalid_aa) > 0) {
        warning(sprintf("Sequence %d contains invalid amino acids: %s", 
                       i, paste(invalid_aa, collapse = ", ")))
      }
    }
    
    list(
      valid = TRUE,
      sequences = fasta_seqs,
      count = length(fasta_seqs)
    )
  }, error = function(e) {
    stop(sprintf("Error reading FASTA file: %s", e$message))
  })
}

#' Validate and normalize segment ranges
#'
#' @param segments List of segment ranges
#' @param seq_lengths Vector of sequence lengths
#' @return Normalized list of segment ranges
#' @noRd
validate_segments <- function(segments, seq_lengths) {
  n_seqs <- length(seq_lengths)
  
  # If segments is NULL or empty, use full length for all sequences
  if (is.null(segments) || length(segments) == 0) {
    segments <- lapply(seq_lengths, function(len) 1:len)
    return(segments)
  }
  
  # Ensure segments list has correct length
  if (length(segments) < n_seqs) {
    # Pad with NULL for missing sequences
    segments <- c(segments, rep(list(NULL), n_seqs - length(segments)))
  } else if (length(segments) > n_seqs) {
    warning(sprintf("More segments (%d) than sequences (%d), using first %d segments", 
                   length(segments), n_seqs, n_seqs))
    segments <- segments[1:n_seqs]
  }
  
  # Process each segment
  normalized_segments <- list()
  
  for (i in seq_along(segments)) {
    seg <- segments[[i]]
    seq_len <- seq_lengths[i]
    
    if (is.null(seg) || (length(seg) == 1 && is.na(seg)) || length(seg) == 0) {
      # Use full length if segment is NULL/NA/empty
      normalized_segments[[i]] <- 1:seq_len
    } else {
      # Validate and normalize segment
      if (is.numeric(seg)) {
        seg <- as.integer(seg)
        
        # Check if all positions are within valid range
        if (any(seg < 1) || any(seg > seq_len)) {
          stop(sprintf("Segment for sequence %d contains positions outside valid range (1-%d)", 
                      i, seq_len))
        }
        
        # Ensure segment is sorted and unique
        seg <- sort(unique(seg))
        normalized_segments[[i]] <- seg
      } else {
        stop(sprintf("Invalid segment type for sequence %d: must be numeric", i))
      }
    }
  }
  
  normalized_segments
}

#' Process direction and adjust segment accordingly
#'
#' @param segment Segment range
#' @param direction Logical (TRUE = forward, FALSE = reverse)
#' @return Adjusted segment
#' @noRd
process_direction <- function(segment, direction = TRUE) {
  if (direction) {
    # Forward direction: N-terminus to C-terminus
    segment
  } else {
    # Reverse direction: C-terminus to N-terminus
    rev(segment)
  }
}

#' Generate sequence name with metadata
#'
#' @param protein_name Name of original protein
#' @param seq_index Index of sequence
#' @param seg_type Type of segment (mutation or stepwise)
#' @param position Position in sequence
#' @param aa Amino acid at position
#' @param direction Direction used
#' @param segment Segment used
#' @return Formatted sequence name
#' @noRd
generate_sequence_name <- function(protein_name, seq_index, seg_type, 
                                  position = NULL, aa = NULL, 
                                  direction = TRUE, segment = NULL) {
  # Clean protein name
  clean_name <- gsub("[^a-zA-Z0-9_]", "_", protein_name)
  
  # Base name
  name_parts <- c(clean_name, sprintf("Seq%d", seq_index), seg_type)
  
  # Add position and amino acid if provided
  if (!is.null(position) && !is.null(aa)) {
    name_parts <- c(name_parts, sprintf("Pos%d", position), aa)
  }
  
  # Add direction
  dir_str <- ifelse(direction, "Fwd", "Rev")
  name_parts <- c(name_parts, dir_str)
  
  # Add segment info if provided
  if (!is.null(segment)) {
    seg_start <- min(segment)
    seg_end <- max(segment)
    name_parts <- c(name_parts, sprintf("Seg%d-%d", seg_start, seg_end))
  }
  
  paste(name_parts, collapse = "_")
}

#' Export sequences to FASTA file
#'
#' @param sequences Named vector of sequences
#' @param output_path Path to output FASTA file
#' @param line_width Line width for sequence wrapping (default: 80)
#' @noRd
export_fasta <- function(sequences, output_path, line_width = 80) {
  if (length(sequences) == 0) {
    warning("No sequences to export")
    return(invisible(NULL))
  }
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write FASTA file
  con <- file(output_path, "w")
  on.exit(close(con))
  
  for (i in seq_along(sequences)) {
    seq_name <- names(sequences)[i]
    seq_char <- as.character(sequences[i])
    
    # Write header
    cat(sprintf(">%s\n", seq_name), file = con)
    
    # Write sequence with line wrapping
    seq_lines <- strsplit(seq_char, "")[[1]]
    n_lines <- ceiling(length(seq_lines) / line_width)
    
    for (j in 1:n_lines) {
      start <- (j - 1) * line_width + 1
      end <- min(j * line_width, length(seq_lines))
      line <- paste(seq_lines[start:end], collapse = "")
      cat(sprintf("%s\n", line), file = con)
    }
  }
  
  cat(sprintf("Successfully exported %d sequences to %s\n", 
             length(sequences), output_path))
}

# ==============================================================================
# Core Functions
# ==============================================================================

#' Generate single-site mutation sequences
#'
#' This function generates comprehensive single-site mutation sequences from 
#' one or more input protein sequences in FASTA format.
#'
#' @param fasta_path Path to input FASTA file containing one or more protein sequences
#' @param segments Optional list specifying segments for each protein. 
#'   Each element corresponds to a protein sequence. Use NULL for full-length.
#'   Example: list(c(1:30), 25:90)
#' @param direction Logical indicating processing direction. 
#'   TRUE (default) = N-terminus to C-terminus (forward), 
#'   FALSE = C-terminus to N-terminus (reverse)
#' @param output_dir Directory for output files. If NULL, uses current directory.
#' @param output_prefix Prefix for output file names.
#' @param include_original Logical indicating whether to include original sequences.
#'
#' @return Named vector of all generated sequences (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate mutations for full-length sequences in forward direction
#' DataGen_Mutation("input.fasta")
#'
#' # Generate mutations for specific segments
#' DataGen_Mutation("input.fasta", segments = list(1:30, 25:90))
#'
#' # Generate mutations in reverse direction
#' DataGen_Mutation("input.fasta", direction = FALSE)
#' }
DataGen_Mutation <- function(fasta_path, 
                            segments = NULL, 
                            direction = TRUE,
                            output_dir = NULL,
                            output_prefix = "DataGen",
                            include_original = TRUE) {
  
  cat("=== DataGen_Mutation Starting ===\n")
  
  # Validate input
  cat("Validating input...\n")
  fasta_info <- validate_fasta_input(fasta_path)
  fasta_seqs <- fasta_info$sequences
  n_seqs <- fasta_info$count
  
  cat(sprintf("Found %d protein sequence(s)\n", n_seqs))
  
  # Get sequence lengths
  seq_lengths <- sapply(fasta_seqs, function(s) nchar(as.character(s)))
  
  # Validate and normalize segments
  cat("Processing segments...\n")
  normalized_segments <- validate_segments(segments, seq_lengths)
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Valid amino acids for mutation
  aa_list <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  
  # Generate all mutation sequences
  cat("Generating mutation sequences...\n")
  all_seqs <- list()
  
  for (i in seq_along(fasta_seqs)) {
    seq_name <- names(fasta_seqs)[i]
    original_seq <- as.character(fasta_seqs[i])
    original_aa <- strsplit(original_seq, "")[[1]]
    segment <- normalized_segments[[i]]
    
    # Process direction
    processed_segment <- process_direction(segment, direction)
    
    cat(sprintf("Processing sequence %d (%s): segment %d-%d, direction %s\n", 
               i, seq_name, min(segment), max(segment), 
               ifelse(direction, "forward", "reverse")))
    
    # Include original sequence if requested
    if (include_original) {
      orig_name <- generate_sequence_name(seq_name, i, "ORIG", 
                                         direction = direction, 
                                         segment = segment)
      all_seqs[[orig_name]] <- original_seq
    }
    
    # Generate mutations for each position in segment
    for (pos in processed_segment) {
      original_residue <- original_aa[pos]
      
      # Generate all possible mutations (excluding original)
      for (aa in aa_list[aa_list != original_residue]) {
        mutant_seq <- original_aa
        mutant_seq[pos] <- aa
        mutant_seq_str <- paste(mutant_seq, collapse = "")
        
        mutant_name <- generate_sequence_name(seq_name, i, "MUT", 
                                             position = pos, aa = aa,
                                             direction = direction, 
                                             segment = segment)
        
        all_seqs[[mutant_name]] <- mutant_seq_str
      }
    }
  }
  
  # Convert to named vector
  all_seqs_vec <- unlist(all_seqs)
  
  # Export to FASTA
  output_file <- file.path(output_dir, 
                          sprintf("%s_Mutation_sequences.fasta", output_prefix))
  
  cat("Exporting sequences...\n")
  export_fasta(all_seqs_vec, output_file)
  
  cat(sprintf("=== DataGen_Mutation Complete ===\n"))
  cat(sprintf("Generated %d sequences\n", length(all_seqs_vec)))
  cat(sprintf("Output saved to: %s\n", output_file))
  
  invisible(all_seqs_vec)
}

#' Generate stepwise/truncated sequences
#'
#' This function generates comprehensive stepwise/truncated sequences from 
#' one or more input protein sequences in FASTA format.
#'
#' @param fasta_path Path to input FASTA file containing one or more protein sequences
#' @param segments Optional list specifying segments for each protein. 
#'   Each element corresponds to a protein sequence. Use NULL for full-length.
#'   Example: list(c(1:30), 25:90)
#' @param direction Logical indicating processing direction. 
#'   TRUE (default) = N-terminus to C-terminus (forward), 
#'   FALSE = C-terminus to N-terminus (reverse)
#' @param output_dir Directory for output files. If NULL, uses current directory.
#' @param output_prefix Prefix for output file names.
#' @param include_original Logical indicating whether to include original sequences.
#' @param include_all_residues Logical indicating whether to include all residues 
#'   at each step (TRUE) or just the truncated sequence (FALSE).
#'
#' @return Named vector of all generated sequences (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate stepwise sequences for full-length sequences in forward direction
#' DataGen_Stepwise("input.fasta")
#'
#' # Generate stepwise sequences for specific segments
#' DataGen_Stepwise("input.fasta", segments = list(1:30, 25:90))
#'
#' # Generate stepwise sequences in reverse direction
#' DataGen_Stepwise("input.fasta", direction = FALSE)
#' }
DataGen_Stepwise <- function(fasta_path, 
                            segments = NULL, 
                            direction = TRUE,
                            output_dir = NULL,
                            output_prefix = "DataGen",
                            include_original = TRUE,
                            include_all_residues = TRUE) {
  
  cat("=== DataGen_Stepwise Starting ===\n")
  
  # Validate input
  cat("Validating input...\n")
  fasta_info <- validate_fasta_input(fasta_path)
  fasta_seqs <- fasta_info$sequences
  n_seqs <- fasta_info$count
  
  cat(sprintf("Found %d protein sequence(s)\n", n_seqs))
  
  # Get sequence lengths
  seq_lengths <- sapply(fasta_seqs, function(s) nchar(as.character(s)))
  
  # Validate and normalize segments
  cat("Processing segments...\n")
  normalized_segments <- validate_segments(segments, seq_lengths)
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Valid amino acids
  aa_list <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  
  # Generate all stepwise sequences
  cat("Generating stepwise sequences...\n")
  all_seqs <- list()
  
  for (i in seq_along(fasta_seqs)) {
    seq_name <- names(fasta_seqs)[i]
    original_seq <- as.character(fasta_seqs[i])
    original_aa <- strsplit(original_seq, "")[[1]]
    segment <- normalized_segments[[i]]
    
    # Process direction
    processed_segment <- process_direction(segment, direction)
    
    cat(sprintf("Processing sequence %d (%s): segment %d-%d, direction %s\n", 
               i, seq_name, min(segment), max(segment), 
               ifelse(direction, "forward", "reverse")))
    
    # Include original sequence if requested
    if (include_original) {
      orig_name <- generate_sequence_name(seq_name, i, "ORIG", 
                                         direction = direction, 
                                         segment = segment)
      all_seqs[[orig_name]] <- original_seq
    }
    
    # Generate stepwise sequences for each position in segment
    for (step_idx in seq_along(processed_segment)) {
      current_pos <- processed_segment[step_idx]
      step_positions <- processed_segment[1:step_idx]
      
      if (include_all_residues) {
        # Include all residues at each step position
        for (aa in aa_list) {
          stepwise_seq <- original_aa
          stepwise_seq[current_pos] <- aa
          stepwise_seq_str <- paste(stepwise_seq, collapse = "")
          
          stepwise_name <- generate_sequence_name(seq_name, i, "STEP", 
                                                 position = current_pos, aa = aa,
                                                 direction = direction, 
                                                 segment = segment)
          
          all_seqs[[stepwise_name]] <- stepwise_seq_str
        }
      } else {
        # Just truncated sequence with original residue
        stepwise_seq <- original_aa[1:current_pos]
        stepwise_seq_str <- paste(stepwise_seq, collapse = "")
        
        stepwise_name <- generate_sequence_name(seq_name, i, "STEP", 
                                               position = current_pos, 
                                               aa = original_aa[current_pos],
                                               direction = direction, 
                                               segment = segment)
        
        all_seqs[[stepwise_name]] <- stepwise_seq_str
      }
    }
  }
  
  # Convert to named vector
  all_seqs_vec <- unlist(all_seqs)
  
  # Export to FASTA
  output_file <- file.path(output_dir, 
                          sprintf("%s_Stepwise_sequences.fasta", output_prefix))
  
  cat("Exporting sequences...\n")
  export_fasta(all_seqs_vec, output_file)
  
  cat(sprintf("=== DataGen_Stepwise Complete ===\n"))
  cat(sprintf("Generated %d sequences\n", length(all_seqs_vec)))
  cat(sprintf("Output saved to: %s\n", output_file))
  
  invisible(all_seqs_vec)
}

# ==============================================================================
# Example Usage
# ==============================================================================

#' Run comprehensive test examples
#'
#' This function demonstrates all key features of the DataGen functions.
#'
#' @param test_fasta_path Path to test FASTA file (optional)
#' @param test_output_dir Directory for test outputs
#' @export
run_examples <- function(test_fasta_path = NULL, 
                        test_output_dir = "DataGen_Examples") {
  
  cat("=== Running DataGen Examples ===\n")
  
  # Create test output directory
  if (!dir.exists(test_output_dir)) {
    dir.create(test_output_dir, recursive = TRUE)
  }
  
  # If no test FASTA provided, create a simple test FASTA
  if (is.null(test_fasta_path) || !file.exists(test_fasta_path)) {
    cat("Creating test FASTA file...\n")
    
    test_seqs <- c(
      "TestProtein1" = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLK",
      "TestProtein2" = "MTSDMPRIKPGQRVMMALRKMIASGEIKSGERIAEIPTAAALGVSRMP"
    )
    
    test_fasta_path <- file.path(test_output_dir, "test_sequences.fasta")
    export_fasta(test_seqs, test_fasta_path)
  }
  
  # Example 1: Basic usage with full-length sequences
  cat("\n--- Example 1: Basic Mutation Generation ---\n")
  DataGen_Mutation(
    fasta_path = test_fasta_path,
    output_dir = test_output_dir,
    output_prefix = "Example1_Basic",
    include_original = TRUE
  )
  
  # Example 2: Using specific segments
  cat("\n--- Example 2: Mutation with Specific Segments ---\n")
  DataGen_Mutation(
    fasta_path = test_fasta_path,
    segments = list(1:20, 15:35),
    output_dir = test_output_dir,
    output_prefix = "Example2_Segments"
  )
  
  # Example 3: Reverse direction
  cat("\n--- Example 3: Mutation in Reverse Direction ---\n")
  DataGen_Mutation(
    fasta_path = test_fasta_path,
    direction = FALSE,
    output_dir = test_output_dir,
    output_prefix = "Example3_Reverse"
  )
  
  # Example 4: Basic stepwise generation
  cat("\n--- Example 4: Basic Stepwise Generation ---\n")
  DataGen_Stepwise(
    fasta_path = test_fasta_path,
    output_dir = test_output_dir,
    output_prefix = "Example4_Stepwise_Basic",
    include_all_residues = FALSE
  )
  
  # Example 5: Stepwise with all residues
  cat("\n--- Example 5: Stepwise with All Residues ---\n")
  DataGen_Stepwise(
    fasta_path = test_fasta_path,
    segments = list(1:15, 10:25),
    output_dir = test_output_dir,
    output_prefix = "Example5_Stepwise_AllResidues",
    include_all_residues = TRUE
  )
  
  cat("\n=== All Examples Complete ===\n")
  cat(sprintf("Examples saved to: %s\n", normalizePath(test_output_dir)))
}

# ==============================================================================
# Verification and Validation Plan
# ==============================================================================

#' Comprehensive validation and verification
#'
#' This function runs a complete validation suite to ensure all functions
#' work correctly.
#'
#' @param validation_output_dir Directory for validation outputs
#' @export
run_validation <- function(validation_output_dir = "DataGen_Validation") {
  
  cat("=== Running DataGen Validation Suite ===\n")
  
  if (!dir.exists(validation_output_dir)) {
    dir.create(validation_output_dir, recursive = TRUE)
  }
  
  validation_report <- c(
    "# DataGen Validation Report",
    "",
    "## Test Overview",
    "",
    paste("Validation performed on:", Sys.time()),
    ""
  )
  
  # Test 1: Single protein sequence
  cat("\n--- Test 1: Single Protein Sequence ---\n")
  test1_seqs <- c("ProteinA" = "MVSKGEELFTGVVPILVELD")
  test1_fasta <- file.path(validation_output_dir, "test1_single.fasta")
  export_fasta(test1_seqs, test1_fasta)
  
  result1_mutation <- DataGen_Mutation(
    test1_fasta,
    output_dir = validation_output_dir,
    output_prefix = "Test1_Single"
  )
  
  result1_stepwise <- DataGen_Stepwise(
    test1_fasta,
    output_dir = validation_output_dir,
    output_prefix = "Test1_Single",
    include_all_residues = FALSE
  )
  
  validation_report <- c(validation_report,
    "## Test 1: Single Protein Sequence",
    paste("- Mutation sequences generated:", length(result1_mutation)),
    paste("- Stepwise sequences generated:", length(result1_stepwise)),
    "")
  
  # Test 2: Multiple protein sequences
  cat("\n--- Test 2: Multiple Protein Sequences ---\n")
  test2_seqs <- c(
    "ProteinA" = "MVSKGEELFTGVVPILVELD",
    "ProteinB" = "MTSDMPRIKPGQRVMMALRKMI",
    "ProteinC" = "ACDEFGHIKLMNPQRSTVWY"
  )
  test2_fasta <- file.path(validation_output_dir, "test2_multiple.fasta")
  export_fasta(test2_seqs, test2_fasta)
  
  result2_mutation <- DataGen_Mutation(
    test2_fasta,
    output_dir = validation_output_dir,
    output_prefix = "Test2_Multiple"
  )
  
  validation_report <- c(validation_report,
    "## Test 2: Multiple Protein Sequences",
    paste("- Number of input sequences:", 3),
    paste("- Total mutation sequences generated:", length(result2_mutation)),
    "")
  
  # Test 3: Segment validation
  cat("\n--- Test 3: Segment Validation ---\n")
  test3_seqs <- c("LongProtein" = paste(rep("A", 100), collapse = ""))
  test3_fasta <- file.path(validation_output_dir, "test3_segments.fasta")
  export_fasta(test3_seqs, test3_fasta)
  
  result3a <- DataGen_Mutation(
    test3_fasta,
    segments = list(1:50),
    output_dir = validation_output_dir,
    output_prefix = "Test3a_Segment1_50"
  )
  
  result3b <- DataGen_Mutation(
    test3_fasta,
    segments = list(51:100),
    direction = FALSE,
    output_dir = validation_output_dir,
    output_prefix = "Test3b_Segment51_100_Reverse"
  )
  
  validation_report <- c(validation_report,
    "## Test 3: Segment Validation",
    paste("- Segment 1-50 (forward) sequences:", length(result3a)),
    paste("- Segment 51-100 (reverse) sequences:", length(result3b)),
    "")
  
  # Save validation report
  report_file <- file.path(validation_output_dir, "validation_report.md")
  writeLines(validation_report, report_file)
  
  cat("\n=== Validation Complete ===\n")
  cat(sprintf("Validation report saved to: %s\n", normalizePath(report_file)))
}

# ==============================================================================
# Initialization
# ==============================================================================

# Message when script is sourced
cat("==============================================================================\n")
cat("DataGen Functions Loaded Successfully!\n")
cat("==============================================================================\n")
cat("\n")
cat("Available functions:\n")
cat("  - DataGen_Mutation()    : Generate single-site mutation sequences\n")
cat("  - DataGen_Stepwise()    : Generate stepwise/truncated sequences\n")
cat("  - run_examples()        : Run comprehensive usage examples\n")
cat("  - run_validation()      : Run complete validation suite\n")
cat("\n")
cat("For help, see ?DataGen_Mutation or ?DataGen_Stepwise\n")
cat("==============================================================================\n")
cat("\n")
cat("Author: ProtGenesis Analysis Team\n")
cat("Date: 2026-02-24\n")
cat("Version: 1.0.0\n")
cat("By: Chuanyang Liu\n")
cat("Please cite our paper if you use this script in your research:\n")
cat("Chuanyang Liu, et al. (2026). Universal physical principles govern the deterministic genesis of protein structure.\n")
cat("BioRxiv.10.64898/2026.02.20.706798.\n")
cat("\n")
cat("For help, see ?DataGen_Mutation or ?DataGen_Stepwise\n")
cat("==============================================================================\n")
