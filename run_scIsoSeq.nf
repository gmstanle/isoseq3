#!/usr/bin/env nextflow

/*
* Input is a single .bam file of CCS reads (name.ccs.bam) plus index file.
* If there is output from multiple flowcells, merge them with bamtools merge 
* (NOT samtools,see https://github.com/PacificBiosciences/PacBioFileFormats/wiki/BAM-recipes#merging).
* Then index with pbindex.
* Future versions of the pipeline should include these steps for reproducibility.
*/

/*
* Most of the processes in the pipeline are designed to save their output, even though mostly their
* output is intermediate and not used for downstream analysis. That is why they have 
* publishDir calls and why everything is saved as a file with a specific name, rather than
* just passed to a channel and therefore only saved to tmp working folders (which would be the most Nextflow way of doing it). 
* This allows for the output of each intermediate step to be analyzed post facto.
* That makes sense since this is an early stage experiment,
* but makes the code less clean and increases the storage required. A pipeline for this protocol,
* once it is more standardized, would avoid the extra code of saving those intermediate files 
* and just pass them directly into channels.
* See https://github.com/nextflow-io/rnatoy for an example pipeline where only the final file
* is saved to an output directory.
*/

// specify on command line:
// params.input (--input <input_dir_path>)
// params.output (--output <output_dir_path>)

//params.merge = true
params.align = true

// Added by Geoff

params.ref_fasta = params.genome ? params.genomes[ params.genome ].fasta_file ?: false : false
params.intron_max = params.genome ? params.genomes[ params.genome ].intron_max ?: false : false
params.primers = params.primer_type ? params.primers_stets[ params.primer_type ].primer_file ?: false : false

log.info "IsoSeq3 NF  ~  version 3.1"
log.info "====================================="
log.info "input paths: ${params.input}"
log.info "output paths: ${params.output}"
log.info "primer set: ${params.primer_type}"
// log.info "merge smrt cells: ${params.merge}"
log.info "align reads: ${params.align}"
log.info "genome: ${params.genome}"
log.info "genome sequence: ${params.ref_fasta}"
log.info "intron max length: ${params.intron_max}"
log.info "\n"


// Geoff: this will need to be changed for ccs input. For now the easiest to do 
//        is just comment out the old channel and instead define the ccs_out channel
//        using fromPath input channels
//Channel
//    // get a pair of bam and bam.pbi files, change file.name to be the 'base' name (w/o .bam or .pbi)
//    // see https://pbbam.readthedocs.io/en/latest/tools/pbindex.html for info on making .pbi file
//    .fromFilePairs(params.input + '*.{bam,bam.pbi}') { file -> file.name.replaceAll(/.bam|.pbi$/,'') }
//    .ifEmpty { error "Cannot find matching bam and pbi files: $params.input. Make sure your bam files are pb indexed." }
//    .set { ccs_in }
//// see https://github.com/nextflow-io/patterns/blob/926d8bdf1080c05de406499fb3b5a0b1ce716fcb/process-per-file-pairs/main2.nf

// Geoff: for using indexed ccs reads as input
Channel:
    .fromFilePairs(params.input + '*.ccs.{bam,bam.pbi}') { file -> file.name.replaceAll(/.ccs.bam$|.ccs.bam.pbi$/,'') }
    .ifEmpty { error "Cannot find matching bam and pbi files: $params.input. Make sure your bam files are pb indexed." }
    .set(ccs_out_indexed)
Channel
    .fromPath(params.input + '*.bam')
    .ifEmpty { error "Cannot find matching bam files: $params.input." }
    .tap { bam_files }

    // make a matching filename 'base' for every file
    .map{ file -> tuple(file.name.replaceAll(/.bam$/,''), file) } 
    .tap { bam_names }

Channel
    .fromPath(params.primers)
    .ifEmpty { error "Cannot find primer file: $params.primers" }
    .into { primers_remove; primers_refine } // puts the primer files into these two channels

// Question: why is the reference a channel? Can't it just be  aregulr 
Channel
    .fromPath(params.ref_fasta)

    // I assume this is a mistake. This should say Cannot find reference file $params.ref_fasta
//    .ifEmpty { error "Cannot find primer file: $params.primers" }
    .ifEmpty { error "Cannot find reference file: $params.ref_fasta" }
    .set {ref_fasta}

// Geoff: why is this step needed? The 
//TODO replace with specific stating of the pbi
Channel
    .fromPath(params.input + '*.bam.pbi')
    .ifEmpty { error "Cannot find matching bam.pbi files: $params.input." }
    .into { pbi_merge_trans; pbi_merge_sub; pbi_polish }

// This process is currently broken since I am using indexed, ccs reads as input.
process ccs_calling{

        tag "circular consensus sequence calling: $name"
        
        publishDir "$params.output/$name/ccs", mode: 'copy'

        input:
        set name, file(bam) from ccs_in.dump(tag: 'input')

        // To make this compatible with the new pipeline, need to add an indexing step
        // and a channel ccs_out_indexed, which is now the input into remove_primers
        // That output will need to be in the same format as .fromFilePairs output
        output:
        file "*"
        set val(name), file("${name}.ccs.bam") into ccs_out
        
        // Geoff
        when:
            !params.input_is_ccs
     
        //TODO make minPasses param as parameter
        """
        ccs ${name}.bam ${name}.ccs.bam --noPolish --minPasses 1
        """
}


// Geoff: TODO: edit this for the primer design I used 
// Geoff: changed to use indexed .bam files (.pbi)
// Geoff: renamed demux which makes for sense for multiplexed samples
//        lima --ccs both demuxes and removes primers
//process remove_primers{
process demux{

    tag "primer removal: $name"

    publishDir "$params.output/$name/lima", mode: 'copy'

    input:
    // weird usage of dump - it is normally for debugging.
//    set name, file(bam) from ccs_out.dump(tag: 'ccs_name')
    set name, file(bam) from ccs_out_indexed
    path primers from primers_remove.collect()
    
    output:
    path "*"
    //set val(name), file("${name}.fl.primer_5p--primer_3p.bam") into primers_removed_out
    // TODO: get file output name
    set val(name), file("${name}.trimmed.bam") into trimmed_out 
 
//    """
//    lima $bam $primers ${name}.fl.bam --isoseq --no-pbi
//    """
    """
    lima --ccs $bam $primers ${name}.trimmed.bam
    """
}

process run_refine{

    tag "refining : $name"
    publishDir "$params.output/$name/refine", mode: 'copy'

    input:
    set name, file(bam) from trimmed_out.dump(tag: 'trimmed')
    path primers from primers_refine.collect()
    

    // flnc = full-length non-concatemer
    output:
    path "*"
    set val(name), file("${name}.flnc.bam") into refine_out
 
    //TODO update input & output channels
    """
    isoseq3 refine $bam $primers ${name}.flnc.bam --require-polya
    """

}


// I am not sure whether the cluster and polish steps are necessary. The PacBio IsoSeq3
// page has them included but Liz Tseng's "best practices for single-cell IsoSeq" does not.
// As of v3.2, both clustering and polishing are performed by IsoSeq cluster
process cluster_reads{

    tag "clustering : $name"
    publishDir "$params.output/$name/cluster", mode: 'copy'

    input:
    set name, file(refined) from refine_out.dump(tag: 'cluster')

    output:
    file "*"
    set val(name), file("${name}.polished.bam") into cluster_out

    """
    isoseq3 cluster ${refined} ${name}.polished.bam
    """
}


// Following best practices here:
// https://github.com/Magdoll/cDNA_Cupcake/wiki/Best-practice-for-aligning-Iso-Seq-to-reference-genome:-minimap2,-deSALT,-GMAP,-STAR,-BLAT

process align_reads{

    tag "mapping : $name"

    // not clear if all files produced by the code or just the files specifid
    // in output are copied
    publishDir "$params.output/$name/minimap2", mode: 'copy'

    input:
   // set name, file(sample) from polish_out.dump(tag: 'align')
    set name, file(clustered_bam) from cluster_out
    path ref from ref_fasta.collect()

    output:    
    path "*.{sorted.sam,log}"
    set name, file("${name}.aln.bam") into aligned_out

    when:
    params.align

    """
    minimap2 $ref $clustered_bam \
        -G $params.intron_max \
        -H \
        -ax splice \
        -C 5 \
        -O6,24 \
        -B4 \
        -uf \
        --secondary=no \
        -t ${task.cpus} > ${name}.aln.sam \
        2> ${name}.log

    """
}


// from https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake:-supporting-scripts-for-Iso-Seq-after-clustering-step#collapse-redundant-isoforms-has-genome
process collapse_isoforms{

    publishDir "$params.output/$name/collapse_isoforms", mode: 'copy'


    input:
        set name, file(aligned_sam) from aligned_out

    output:
        path "*{gff,fq,txt}"


    """
    sort -k 3,3 -k 4,4n $aligned_sam > sorted.sam
    collapse_isoforms_by_sam.py --input sorted.sam \
      -s flnc.fasta.sorted.sam -c 0.99 -i 0.95 -o flnc.5merge

    """
}
