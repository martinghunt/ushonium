# ushonium

usher + taxonium on Covid sequences.

## Build

Build container with:

```
singularity build ushonium.img Singularity.def
```


## Usage

The container has a script called `ushonium`, which makes a taxonium
`jsonl.gz` file from fasta consensus sequences.

It uses mafft to align all sequences to the Covid reference (ignoring
indels, so all the same length), makes an optimized tree with usher, and
then uses taxoniumtools to make the taxonium jsonl file.

The filenames of consensus sequences need to be in a tab-delimited file
with no headings and two columns:

1. Name of sample
2. FASTA file.

Assuming that TSV file is called `samples.tsv`, the usage is:

```
ushonium samples.tsv outdir
```

where the output directory `outdir` will be created. The final taxonium file
is called `05.taxonium.jsonl.gz`.

### Use >1 CPU

Most of the time is spent making the MSA, which by default uses 1 cpu.
Run in `N` in parallel using the option
```
--cpus N
```
This will hugely speed up the script.

### Metadata

You can set the title that will appear in the taxonium browser with:
```
--title "My awesome tree"
```

You can also use a file of metadata for each sample (that has eg lineage,
country etc). This needs to be tab-delimited with the first line having
column headers. One column must have the
name of the sample, and must exactly match the name given in `samples.tsv`.
By default, this column is assumed to have the name `strain`, but you
can change it to eg `my_names` with:
```
--metacol_name my_names
```
If the metadata file is called `metadata.tsv`, and we want columns
`col11` and `col2`, then use:
```
--metadata_tsv metadata.tsv --metacols col1,col2
```

### Reference genome

The script `ushonium`  was set up for covid,
using the recommended genbank reference from taxoniumtools.
This genbank file is included in the container, so there is
no need to specify the reference genome when running `ushonium` unless
you want to use a different reference. The option to change it to `my_ref.gb` is
`--ref_gb my_ref.gb`.

From the taxoniumtools documentation: "Right now Taxoniumtools is limited in
the types of genome annotations it can support, for SARS-CoV-2 we recommend
using the exact modified .gb file we use in the example, which splits ORF1ab
into ORF1a and ORF1b to avoid the need to model ribosome slippage."
See https://docs.taxonium.org/en/latest/taxoniumtools.html.
