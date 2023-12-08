# MisMatchFinder
The MisMatchFinder algorithm identifies mismatches within reads compared to the reference genome and filters background noise unrelated to somatic mutations through:
 * the use of high thresholds for mapping and base quality,
 * strict consensus between overlapping read-pairs,
 * gnomAD-based germline variant filtering​,
 * and a ctDNA-centric fragmentomics filter


# Installation

Install time is 100% dependent on the compile time and therefore different between systems, but on a standard x86 laptop compile time is typically less than 2 minutes.

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

# Example analysis

With the [example bam](example/example.bam) runing the default analysis with the provided [whitelist](example/whitelist.bed) we generate the vcf output [`example_bamsites.vcf.gz`](example/example_bamsites.vcf.gz)

please refer to [echtvar release](https://github.com/brentp/echtvar/releases/tag/v0.1.9) to download the echtvar reference for gnomad v3.1.2 

Runtime of this step should be less than 20 seconds on the test data
```
$ mismatchfinder --germline_file /path/to/echtvarfile/gnomad.v3.1.2.echtvar.v2.zip  --whitelist_bed /path/to/whitelist.bed -o /path/to/outputfolder/--strict_overlaps --only_overlaps /path/to/example.bam 
[...]
2023-12-08T04:15:16.586Z INFO [mismatchfinder] Found 8738 mismatches 
2023-12-08T04:15:18.389Z INFO [mismatchfinder] Found 828 somatic mismatches
```

once we have the vcf (can be found in the [examplefolder](example)), we can perform signature deconvolution in R once both data.table and sigminer are installed

```R
library(data.table)
library(sigminer)

maf <- sigminer::read_vcf("example_bamsites.vcf.gz", genome_build = "hg38")
tally <- sigminer::sig_tally(maf, mode="SBS", ref_genome="BSgenome.Hsapiens.UCSC.hg38")

# get the relative signature weights
SBSfit <- sigminer::sig_fit(t(tally$nmf_matrix), sig_index = signatures_selection, sig_db = "SBS", return_class = "data.table", mode="SBS", type="relative")

#plot the signature weights (highlighting SBS7a)
plot(1:ncol(SBSfit[,-1]), round(SBSfit[1,-1], digits=5), type='h', lwd=4, lend=1, col=ifelse(colnames(SBSfit[1,-1])=="SBS7a", "red", "black"), xlab="", ylab="Signature weight", las=2, xaxt='n', bty='n')
text(1:ncol(SBSfit[,-1]), line2user(0,1), labels = colnames(SBSfit[,-1]), srt=45, adj=1, cex=0.5, xpd=TRUE, col=ifelse(colnames(SBSfit[1,-1])=="SBS7a", "black", "lightgrey"))


```

the output on the test data will look similar to this
![plot](example/example_signature_weights.png)



# Known issues and possible extentions

 * No CRAM support yet
 * Zarr support only rudimentary
