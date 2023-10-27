# MisMatchFinder
The MisMatchFinder algorithm identifies mismatches within reads compared to the reference genome and filters background noise unrelated to somatic mutations through:
 * the use of high thresholds for mapping and base quality,
 * strict consensus between overlapping read-pairs,
 * gnomAD-based germline variant filtering​,
 * and a ctDNA-centric fragmentomics filter


# Installation

Requirements: Rust

1. Clone this repository
2. Build the release optimised binary with ```cargo build -r```
3. Move the binary (```target/release/mismatchfinder```) into the ```$PATH```

# Usage

To ensure optimal results, please use the gnomAD echtvar (available from [https://github.com/brentp/echtvar/releases](https://github.com/brentp/echtvar/releases)) or an eqivalent echtvar file

```
USAGE:
    mismatchfinder [OPTIONS] --output <OUTPUT_FOLDER> [--] [BAMS]...

ARGS:
    <BAMS>...    Bams to analyse

OPTIONS:
        --blacklist_bed <BLACKLIST_FILE>
            Bed file for genomic regions to ignore in the analysis. White list bed file regions
            overwrite black list regions

        --fragment_length_intervals <FRAGMENT_LENGTH_INTERVALS>...
            Length of fragments to be considered in the analysis [default: 100-150 250-325]

        --germline_file <GERMLINE_FILE>
            File to read germline information from default is an echtvar file

    -h, --help
            Print help information

        --maximum_edit_distance_per_read <MAX_EDIT_DISTANCE_PER_READ>
            Maximum mismatches we allow in the read [default: 15]

        --minimum_average_base_quality <MIN_AVERAGE_BASE_QUALITY>
            Mimimum average base quality of the read [default: 25]

        --minimum_edit_distance_per_read <MIN_EDIT_DISTANCE_PER_READ>
            Mimimum mismatches we require in the read [default: 1]

    -o, --output <OUTPUT_FOLDER>
            Output folder to write files to

        --only_overlaps
            only use the overlap of the two reads

        --overwrite
            overwrite previous results

    -q, --minimum_mapping_quality <MIN_MAPPING_QUALITY>
            Mimimum mapping quality for a read to be considered [default: 20]

    -Q, --minimum_base_quality <MIN_BASE_QUALITY>
            Mimimum base quality of the mismatch (BQ is summed for a readpair if reads agree)
            [default: 65]

        --strict_overlaps
            only analyse mismatches if the read pair agree (does not restrict to only overlap)

    -V, --version
            Print version information

        --whitelist_bed <WHITELIST_FILE>
            Bed file for genomic regions to include in the analysis. White list bed file regions
            overwrite black list regions [default: all regions]
```

# Known issues and possible extentions

 * No CRAM support yet
 * Zarr support only rudimentary
