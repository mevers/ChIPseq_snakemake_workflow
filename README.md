# General workflow

It is advised not to mix ChIP-seq experiments using different anybodies in the same `snakemake` workflow. In other words, every ChIP-seq experiment with a specific antibody should get its own `snakemake` workflow.

# Naming requirements

## In `config.yaml`

It is important that values of `samples` matchkeys of `units`; for example, consider the following `units` entries

```yaml
units:
    control_WT_rep1:
    	- reads_R1.fastq.gz
        - reads_R2.fastq.gz
```

Then there needs to be a matching value `control_WT_rep1` in `samples`. e.g.

```yaml
samples:
    WT:
        control:
            - control_WT_rep1
```

Notice also that the keys `WT` and `control` have to match `samplea` key.

The general format is

```yaml
samples:
    <CONDITION1>:
        control:
            - control_<CONDITION1>_<REPLICATE>
            - ...
        IP:
            - IP_<CONDITION1>_<REPLICATE>
            - ...
    <CONDITION2>:
        control:
            - control_<CONDITION2>_<REPLICATE>
            - ...
        IP:
            - IP_<CONDITION2>_<REPLICATE>
            - ...
    ...
```

# How to get started

## Changes to `config.yaml` to reflect the actual experiment design and input files

### Samples

Change `samples` and `units` (see the previous section for naming requirements)

### `deeptools` parameters

- Confirm `ignore_for_normalisation` entries
- Set `binsize` for BigWig plots
- Specify `scale_method` which must be one of `["readCount", "SES", "None"]`
- Adjust `bamCompare_specific_args`; e.g. arguments pertaining to the scaling should go here, e.g.

    ```yaml
    deeptools:
        ...
        scale_method: "SES"
        bamCompare_specific_args:
            - "--sampleLength 2000"
        ...
    ```

**A note on scaling**  
Which scaling method is best depends on the data and requirements. For example, if coverage is low `scale_method: "SES"` (which corresponds to `--scaleFactorsMethod SES` in `bamCompare`) might not work, or might require a larger `bamCompare` `--sampleLength` value. In that case, the `readCount` scaling method may be used

```yaml
deeptools:
    ...
    scale_method: "readCount"
    bamCompare_specific_args:
    ...
```

### `macs2` parameters

- Adjust the `format` of the input files; for paired-end data

    ```yaml
    macs2:
        ...
        format: "BAMPE"
    ```
    For single end data

    ```yaml
    macs2:
        ...
        format: "BAM"
    ```

    From the [MACS2 Usage](https://github.com/taoliu/MACS#usage)

    > Format of tag file can be ELAND, BED, ELANDMULTI, ELANDEXPORT, ELANDMULTIPET (for pair-end tags), SAM, BAM, BOWTIE, BAMPE or BEDPE. Default is AUTO which will allow MACS to decide the format automatically. AUTO is also useful when you combine different formats of files. Note that MACS can't detect BAMPE or BEDPE format with AUTO, and you have to implicitly specify the format for BAMPE and BEDPE.
