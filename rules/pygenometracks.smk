import re

def get_pygenometracks_ini(wildcards):
    if re.search(r'U13369.1', wildcards.region):
        return expand(join(
            config["analysisdir"],
            "{reference_version}/pygenometracks",
            "tracks_U13369.1_{operation}_bw{binwidth}.ini"),
            reference_version = wildcards.reference_version,
            operation = wildcards.operation,
            binwidth = wildcards.binwidth)
    else:
        return expand(join(
            config["analysisdir"],
            "{reference_version}/pygenometracks",
            "tracks_genome_{operation}_bw{binwidth}.ini"),
            reference_version = wildcards.reference_version,
            operation = wildcards.operation,
            binwidth = wildcards.binwidth)



rule tidy_IDRsummarised_narrowPeak:
    version: "1.0"
    input:
        join(
            config["analysisdir"],
            "{reference_version}/idr",
            "{treatment}_IDRsummarised.narrowPeak")
    output:
        join(
            config["analysisdir"],
            "{reference_version}/pygenometracks",
            "{treatment}_IDRpeaks.narrowPeak")
    log:
        join(
            "logs",
            "pygenometracks_tidy_IDRsummarised_narrowPeak_" +
            "{reference_version}_{treatment}.log")
    params:
        maxWidth = 10000
    script:
        "../scripts/tidy_idr_narrowPeak.py"


rule create_pygenometracks_ini:
    version: "1.0"
    input:
        bw = expand(join(
            config["analysisdir"],
            "{{reference_version}}/deeptools/bamCompare/{{operation}}",
            "{name}_vs_pooled_control.normSES.bw{{binwidth}}.bw"),
            name = [sample
                for treatment in config["samples"].keys()
                    for sample in config["samples"][treatment]["IP"]]),
        bed_peaks = expand(join(
            config["analysisdir"],
            "{{reference_version}}/pygenometracks",
            "{treatment}_IDRpeaks.narrowPeak"),
            treatment = config["samples"].keys()),
        bed_annot = join(
            config["refdir"],
            config["reference"]["id"],
            config["reference"]["bed"])
    output:
        join(
            config["analysisdir"],
            "{reference_version}/pygenometracks",
            "tracks_{region}_{operation}_bw{binwidth}.ini")
    log:
        join(
            "logs",
            "make_pygenometracks_ini_{reference_version}_{region}_" +
            "{operation}_bw{binwidth}.log")
    params:
        range = [-5, 5]
    script:
        "../scripts/make_pygenometracks_ini.py"


rule pygenometracks_plot_ratio:
    version: "1.0"
    input:
        ini = get_pygenometracks_ini
    output:
        join(
            config["analysisdir"],
            "{reference_version}/pygenometracks",
            "{operation}_IP_vs_pooled_control_{region}_bw{binwidth}.{suffix}")
    log:
        join(
            "logs",
            "pygenometracks_{reference_version}_{operation}_" +
            "IP_vs_pooled_control_{region}_bw{binwidth}_{suffix}.log")
    params:
        cmd = "pyGenomeTracks"
    shell:
        """
            {params.cmd} \
            --tracks {input.ini} \
            --region {wildcards.region} \
            --width 40 \
            --dpi 300 \
            --fontSize 8 \
            --outFileName {output} > {log} 2>&1
        """
