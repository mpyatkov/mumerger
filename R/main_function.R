#' Run the full mumerge peak merging workflow
#'
#' This is the main wrapper function that takes file paths as input,
#' reads the data, processes it in parallel, and returns the final set
#' of merged peaks.
#'
#' @param files_path Path to the directory containing BED/XLS peak files.
#' @param ncores Number of cores to use for parallel processing.
#' @param chunk_size Number of regions to process in each parallel chunk.
#'
#' @return A tibble with the final merged peaks (seqnames, start, end).
#' @export
run_peak_merging <- function(files_path, pattern = "bed", ncores = 1, chunk_size = 500) {

  files_list <- list.files(path = files_path, pattern = pattern, recursive = T, full.names = T)
  extension <- tools::file_ext(files_list[[1]])

  combined_data_granges <- if (extension == "xls") {
    vroom::vroom(files_list, id = "filepath", col_names = T, show_col_types = F, comment = "#") %>%
      dplyr::select(seqnames = chr, start, end)
  } else {
    vroom::vroom(files_list, id = "filepath", col_names = F, show_col_types = F, comment = "#") %>%
      dplyr::select(seqnames = X1, start = X2, end = X3)
  }

  combined_data_granges <- GenomicRanges::makeGRangesFromDataFrame(combined_data_granges, keep.extra.columns = T)

  union <- combined_data_granges %>%
    GenomicRanges::reduce(., min.gapwidth = 1L)

  overlap <- GenomicRanges::findOverlaps(union, combined_data_granges)
  union_overlapping <- union[queryHits(overlap)]
  test_overlapping <- combined_data_granges[subjectHits(overlap)]

  result_table <- dplyr::tibble(
    seqnames = as.character(seqnames(union_overlapping)),
    start_union = start(union_overlapping),
    end_union = end(union_overlapping),
    start = start(test_overlapping),
    end = end(test_overlapping),
    sample_id = test_overlapping$sample_id,
    group_id = test_overlapping$group_id
  )

 singletones <- result_table %>%
    dplyr::add_count(seqnames,start_union,end_union) %>%
    dplyr::filter(n == 1) %>%
    dplyr::select(seqnames,start,end) %>%
    dplyr::arrange(seqnames,start)

  others <- result_table %>%
    dplyr::add_count(seqnames,start_union,end_union) %>%
    dplyr::filter(n > 1) %>%
    dplyr::arrange(seqnames,start)
  #others <- others %>% slice_head(n = 10000)

  if(ncores == 1) {
    chunk_size <- 0

    final_bed_results <- others %>%
      dplyr::mutate(region_id = paste(seqnames, start_union, end_union, sep = "_")) %>%
      dplyr::group_by(region_id) %>%
      tidyr::nest() %>%
      # Using future_map for parallel, but map works for testing.
      dplyr::mutate(new_peaks = map(data, ~process_one_union_region(.x))) %>%
      dplyr::ungroup() %>%
      dplyr::select(new_peaks) %>%
      tidyr::unnest(cols = new_peaks)
  } else {
    # 1. Create the region_id for grouping, but DON'T group the main table yet.
    others_with_id <- others %>%
      dplyr::mutate(region_id = paste(seqnames, start_union, end_union, sep = "_"))

    # 2. Get a list of the unique units of work (the regions).
    unique_region_ids <- unique(others_with_id$region_id)

    # 3. Split the LIST of unique region IDs into chunks for the workers.
    #    This is the key change. We are splitting the groups, not the raw data rows.
    #    Let's set a chunk size for the number of regions per worker, e.g., 1000 regions per chunk.
    regions_per_chunk <- chunk_size     # You can tune this
    chunked_region_ids <- split(unique_region_ids, ceiling(seq_along(unique_region_ids) / regions_per_chunk))

    # 4. Define the function that a worker will execute.
    #    It takes a list of region_ids, filters the main data for those, and then processes them.
    process_region_chunk <- function(region_ids_chunk, full_data) {
      full_data %>%
        # Filter the big table to get only the data for the regions in this chunk
        dplyr::filter(region_id %in% region_ids_chunk) %>%
        # Now, perform the same grouping and processing as the sequential version
        dplyr::group_by(region_id) %>%
        tidyr::nest() %>%
        dplyr::mutate(new_peaks = map(data, ~process_one_union_region(.x))) %>%
        dplyr::ungroup() %>%
        dplyr::select(new_peaks) %>%
        tidyr::unnest(cols = new_peaks)
    }

    # 5. Execute in parallel using future_map
    future::plan(multicore, workers = ncores)

    # We pass the chunks of IDs to the workers.
    # The .options argument is crucial to make the large 'others_with_id' dataframe
    # available to each worker without repeatedly sending it.
    res <- furrr::future_map(
      chunked_region_ids,
      ~process_region_chunk(.x, full_data = others_with_id),
      .options = furrr_options(
        packages = c("tidyverse", "GenomicRanges"),
        globals = c("others_with_id", "process_one_union_region", "process_one_condition",
                    "find_local_maxima_vectorized", "collision_resolver")
      )
    )

    future::plan(sequential)
    final_bed_results <- bind_rows(res)
  }


  final_res <- rbind(singletones, final_bed_results) %>%
    dplyr::arrange(seqnames, start, end)

  return(final_res)
}
