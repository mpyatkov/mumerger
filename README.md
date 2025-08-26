# muMergeR: A Statistical Peak Merging Algorithm in R

Originaly presented here [Transcription factor enrichment analysis (TFEA) quantifies the activity of multiple transcription factors from a single experiment](https://www.nature.com/articles/s42003-021-02153-7), link to original github repo: [https://github.com/Dowell-Lab/TFEA](https://github.com/Dowell-Lab/TFEA)

muMergeR is an R implementation of the muMerge algorithm for intelligently merging and refining genomic peak calls from multiple experiments (e.g., ATAC-seq, ChIP-seq).

Instead of simply taking the union of intervals, this method models each peak as a Gaussian distribution. It then combines the signals from all samples across different experimental conditions into a unified probability landscape. New, consensus peak locations are identified from the local maxima of this landscape, resulting in a robust master peak list that statistically resolves conflicts in overlapping regions.

### Algorithm description

1.  **Model**: Each individual peak (`start`, `end`) is converted into a Gaussian (normal) distribution with a mean (`mu`) and standard deviation (`sigma`).
2.  **Combine**: For simplicity, consider a single isolated region that contains multiple peaks from different conditions. Within this region, the peak probability distributions are summed to create a combined probability landscape. The key difference from the original algorithm is that all information about condition is ignored, and each peak is assigned an equal, condition-independent weight. This resolve the issue when number of replicates are different for each condition [issue](https://github.com/Dowell-Lab/TFEA/issues/33). 
3.  **Find Maxima**: The algorithm identifies the locations of local maxima in the combined landscape, which represent the new consensus peak centers.
4.  **Refine & Resolve**: A new width (`sigma`) is calculated for each new peak based on a weighted average of nearby original peaks. Any remaining overlaps between the new peaks are resolved to produce the final, non-overlapping master peak list.

### Algorithm Visualization

![Illustration of the muMerge algorithm](./docs/comparison.jpg)

> *Figure: Comparison original mumerge implemented in TFEA and current implementation.*

### Perfomance optimization

This implementation supports multithreading, so feel free to use the _ncores_ parameter of the _mumerge::run_peak_merging_ function to speed up the results.


---

## Usage Example in R

This example demonstrates the full workflow, from loading input files to generating the final consensus peak BED file.

```r
install.packages("remotes")   # if not already installed
remotes::install_github("mpyatkov/mumerger", dependencies = TRUE)

# Load necessary libraries
suppressPackageStartupMessages(library(mumerge))
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(tidyverse)) 
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(future))

# --- Argument Parsing ---
p <- arg_parser('Postprocessing of ChIP-seq peaks with mumerge')
p <- add_argument(p, '--files_path', help="Path to MACS2 bed/xls files.")
p <- add_argument(p, '--pattern', default = "narrow", help="File name pattern")
p <- add_argument(p, '--output_file', default="mumerge_output.bed", help="Output BED file.")
p <- add_argument(p, '--ncores', default=1, help="Number of cores for processing.")
p <- add_argument(p, '--chunk_size', default=500, help="Number of regions per chunk for parallel processing.")
argv <- parse_args(p)

# --- Main Logic ---
final_results <- mumerge::run_peak_merging(
  files_path = argv$files_path,
  pattern = argv$pattern,
  ncores = as.integer(argv$ncores),
  chunk_size = as.integer(argv$chunk_size)
)

# --- Write Output ---
vroom::vroom_write(final_results, argv$output_file, col_names = FALSE)

cat("Mumerge processing complete. Results saved to:", argv$output_file, "\n")
```
