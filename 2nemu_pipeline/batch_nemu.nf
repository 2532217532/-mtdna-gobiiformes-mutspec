// batch.nf - 批量运行 NEMU

nextflow.enable.dsl = 2

params.input_dir = "/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/nemu_input"
params.output_dir = "/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/nemu_output"
params.db_base = "/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/ncbi_db"
params.taxdump = "$HOME/.taxonkit"
params.gencode = 2

workflow {
    
    all_fastas = channel.fromPath("${params.input_dir}/*/*.fasta")
        .map { fasta_file ->
            def gene = fasta_file.parent.name
            def base = fasta_file.baseName
            def species = base.replaceAll('.*__', '').replaceAll('_', ' ')
            def out_dir = "${params.output_dir}/${gene}/${base}"
            def db_path = "${params.db_base}/${gene}_BLAST"
            
            [base, fasta_file, species, gene, out_dir, db_path]
        }
    
    pending = all_fastas.filter { base, fasta, species, gene, out_dir, db_path ->
        !file("${out_dir}/ms12syn_labeled.txt").exists()
    }
    
    RUN_NEMU(pending)
}

process RUN_NEMU {
    tag { "$gene/$base" }
    cpus 4
    errorStrategy 'ignore'  // ✅ 改这里
    
    input:
    tuple val(base), path(fasta), val(species), val(gene), val(out_dir), val(db_path)
    
    output:
    tuple val(gene), val(base), optional: true
    
    script:
    """
    mkdir -p ${out_dir}
    
    echo "🚀 处理: $gene/$base"
    
    nextflow run /home/zengjl/bio_soft/nemu-pipeline-nf/main.nf \
        -c /home/zengjl/bio_soft/nemu-pipeline-nf/nextflow.config \
        -process.cpus=4 \
        -w ${out_dir}/work \
        --resume \
        -with-trace ${out_dir}/trace.txt \
        -output-dir ${out_dir} \
        --input_type protein \
        --input $fasta \
        --species_name "$species" \
        --gencode 2 \
        --db $db_path \
        --taxdump ${params.taxdump} \
        > ${out_dir}/nemu.log 2>&1
    """
}