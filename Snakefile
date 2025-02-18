from pathlib import Path
import re
import polars as pl
import pyfastx

# COMMON ------------------------------------------------------------------------------------------------------------------------------------------------- COMMON
def get_fastq(wildcards):
    return [str(file) for file in Path(SAMPLE_DIR_FASTQ).iterdir() if re.search(r"fastq|fq", file.name) and wildcards.sample in file.name]

def get_sample_names_fastq():
    return [file.stem.split(".")[0] for file in Path(SAMPLE_DIR_FASTQ).iterdir() if re.search(r"fastq|fq", file.name)]

# IO ------------------------------------------------------------------------------------------------------------------------------------------------- IO

configfile: "config_16s_vsearch.yaml"


SAMPLE_DIR_FASTQ = config["SAMPLE_DIR_FASTQ"]
SAMPLES_FASTQ = get_sample_names_fastq()
RESULTS = f"{config['RESULTS']}/{Path(SAMPLE_DIR_FASTQ).resolve().parts[-1]}"


THREADS = config["THREADS"]

# RULE ALL ------------------------------------------------------------------------------------------------------------------------------------------------- RULE ALL

rule all:
    input:
        expand(f"{RESULTS}/{{sample}}/usearch/{{sample}}.with_sequences.csv", sample=SAMPLES_FASTQ),
        #expand(f"{RESULTS}/{{sample}}/usearch/{{sample}}.rel_abund.csv", sample=SAMPLES_FASTQ),
        expand(f"{RESULTS}/{{sample}}/usearch/{{sample}}.rel_abund_from_all_reads.csv", sample=SAMPLES_FASTQ),

# RULES ------------------------------------------------------------------------------------------------------------------------------------------------- RULES

#rule chopper:
#    input:


rule seqkit_rmdup:
    input:
        fastq=get_fastq
    output:
        deduped=f"{RESULTS}/{{sample}}/deduped/{{sample}}.deduped.fastq"
    shell:
        """
        zcat {input.fastq} | seqkit rmdup --by-seq -o {output.deduped}
        """

rule seqkit_fq2fa:
    input:
        fastq=rules.seqkit_rmdup.output.deduped,
    output:
        fasta=f"{RESULTS}/{{sample}}/fasta/{{sample}}.fasta",
    shell:
        """
        cat {input.fastq} | seqkit fq2fa -o {output.fasta}
        """
        
rule vsearch_chimera:
    input:
        fasta=rules.seqkit_fq2fa.output.fasta,
    params:
        gold=config["CHIMERA"],
        threads=config["THREADS"],
    log:
        chimera_log=f"{RESULTS}/{{sample}}/logs/{{sample}}.chimera.log",
    output:
        no_chimeras=f"{RESULTS}/{{sample}}/no_chimeras/{{sample}}.no_chimeras.fasta",
    shell:
        """
        vsearch --uchime_ref {input.fasta} --db {params.gold} --nonchimeras {output.no_chimeras} --threads {params.threads} > {log.chimera_log} 2>&1
        """
    
    
rule vsearch_usearch:
    input:
        fasta=rules.vsearch_chimera.output.no_chimeras,
    output:
        blast6=f"{RESULTS}/{{sample}}/usearch/{{sample}}.blast6",
        uc=f"{RESULTS}/{{sample}}/usearch/{{sample}}.uc",
    log:
        usearch_log=f"{RESULTS}/{{sample}}/logs/{{sample}}.usearch.log",
    params:
        gtdb=config["GTDB"],
        threads=config["THREADS"],
        id_percent=config["ID_PERCENT"],
        max_accepts=config["MAX_ACCEPTS"],
    shell:
        """
        vsearch --usearch_global {input.fasta} --db {params.gtdb} --id {params.id_percent} --blast6out {output.blast6} --uc {output.uc} --threads {params.threads} --strand both --maxaccepts {params.max_accepts} > {log.usearch_log} 2>&1
        """

rule parse_vsearch_blast6:
    input:
        blast6=rules.vsearch_usearch.output.blast6,
    output:
        parsed_blast6=f"{RESULTS}/{{sample}}/usearch/{{sample}}.blast6.csv",
    log:
        parse_blast6=f"{RESULTS}/{{sample}}/logs/{{sample}}.parse_blast6.log",
    run:
        def parse_blast6(file: Path):
            columns = [
                "seq_name",
                "species",
                "identity",
                "alignment_len",
                "mismatches",
                "gaps",
                "start_query",
                "end_query",
                "start_target",
                "end_target",
                "e_value",
                "bit_score"
            ]

            raw = (
                pl.read_csv(file, separator="\t", has_header=False, new_columns=columns)
                # remove all numbers from species name
                .filter(~pl.col("identity").is_null())
                .with_columns(
                    species=pl.when(pl.col("species").str.contains(r"\d")).then(pl.col("species").str.replace_all(r"\d", "")).otherwise(pl.col("species")) 
                )
                .with_columns(species=pl.col("species").str.replace_all(r"_\w{1}_", "_"))
                .with_columns(species=pl.col("species").str.replace_all(r"_\w{1}$", ""))
                .with_columns(count=pl.col("species").n_unique().over("seq_name"))
                .filter(pl.col("species") != "None")
            )

            total_sequences = raw.n_unique("seq_name")
            one_sequences = raw.filter(pl.col("count") == 1).n_unique("seq_name")
            two_sequences = raw.filter(pl.col("count") == 2).n_unique("seq_name")
            three_sequences = raw.filter(pl.col("count") == 3).n_unique("seq_name")

            with open(log.parse_blast6, 'a+') as log_file:
                print(f"Total sequences: {total_sequences}", file=log_file)
                print(f"Sequences classified as 1 species: {one_sequences}", file=log_file)
                print(f"Sequences classified as 2 species: {two_sequences}", file=log_file)
                print(f"Sequences classified as 3 species: {three_sequences}", file=log_file)

            df = (
                raw
                .sort(["seq_name", "identity", "mismatches"], descending=[True, True, False])
                .unique("seq_name", keep="first", maintain_order=True)
                .unique("seq_name")
                .drop(["bit_score", "e_value", "start_target", "end_target"])
                .with_columns(file=pl.lit(file.stem))
            )
            return df
        
        df = parse_blast6(Path(input.blast6))
        df.write_csv(output.parsed_blast6)

rule make_relative_abundance:
    input:
        blast6=rules.parse_vsearch_blast6.output.parsed_blast6
    output:
        rel_abund=f"{RESULTS}/{{sample}}/usearch/{{sample}}.rel_abund_from_all_reads.csv",
    params:
        count_filter=config["COUNT_FILTER"],
    run:
        def make_rel_abund(file, count_filter):
            df = (
                pl.read_csv(file)
                .with_columns(total_reads=pl.len())
                .filter(pl.col("count") <= count_filter)
                .group_by(["species", "file", "total_reads"])
                .agg(
                    n_reads=pl.len()    
                )
                .sort("n_reads", descending=True)
                .with_columns(rel_abund=pl.col("n_reads") / pl.col("total_reads"))
                ####
            )

            return df
        
        rel_abund = make_rel_abund(input.blast6, params.count_filter)
        rel_abund.write_csv(output.rel_abund)


rule add_sequence:
    input:
        blast6=rules.parse_vsearch_blast6.output.parsed_blast6,
        fasta=rules.vsearch_chimera.output.no_chimeras,
    output:
        added_sequences=f"{RESULTS}/{{sample}}/usearch/{{sample}}.with_sequences.csv",
    run:
        def read_fasta(fasta):
            fasta = pyfastx.Fastx(fasta)
            name = []
            seq = []
            for x in fasta:
                name.append(x[0])
                seq.append(x[1])

            return pl.DataFrame({"seq_name": name, "seq": seq})


        fasta_df = read_fasta(input.fasta)
        df = pl.read_csv(input.blast6)
        
        joined = df.join(fasta_df, on="seq_name")
        joined.write_csv(output.added_sequences)

    
