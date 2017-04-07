import "mutect2.wdl" as m2

workflow create_m2_pon {
    File normal_bam_list
    File normal_bai_list
    Array[File] normal_bams = read_lines(normal_bam_list)
    Array[File] normal_bais = read_lines(normal_bai_list)
    Array[Pair[File, File]] bam_pairs = zip(normal_bams, normal_bais)
    File intervals
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File cosmic
    File cosmic_index
    File picard_jar
    String m2_docker
    File? gatk4_jar_override
    String output_name

    scatter (p in bam_pairs) {
         call m2.Mutect2 as m2_full {
             input:
                 gatk4_jar="/root/gatk-protected.jar",
                 intervals = intervals,
                 ref_fasta = ref_fasta,
                 ref_fasta_index = ref_fasta_index,
                 ref_dict = ref_dict,
                 tumor_bam = p.left,
                 tumor_bam_index = p.right,
                 tumor_sample_name = sub(sub(p.left, "[gs://]*[/]*.*/", ""), "\\.bam$", ""),
                 pon_index = pon_index,
                 scatter_count = 50,
                 cosmic = cosmic,
                 cosmic_index = cosmic_index,
                 is_run_orientation_bias_filter = false,
                 is_run_oncotator = false,
                 oncotator_docker = "dummy",
                 m2_docker = m2_docker,
                 gatk4_jar_override = gatk4_jar_override,
                 preemptible_attempts = 2,
                 artifact_modes = [],
                 picard_jar = picard_jar
         }
    }

    call CreateSomaticPanelOfNormals {
        input:
             normal_vcfs=m2_full.unfiltered_vcf,
             output_name=output_name
    }

    output {
        File m2_pon = CreateSomaticPanelOfNormals.m2_pon
    }
}

task CreateSomaticPanelOfNormals {
    Array[File] normal_vcfs
    String output_name

    command {
        # java -jar gatk-protected.jar CreateSomaticPanelOfNormals -vcfs foo.vcf -vcfs foo2.vcf -O foo.pon.vcf
        java -jar gatk-protected.jar CreateSomaticPanelOfNormals -vcfs ${sep=" -vcfs " normal_vcfs} -O ${output_name}.pon.vcf
    }
    runtime {
        docker: "${m2_docker}"
        memory: "12 GB"
        disks: "local-disk " + 500 + " SSD"
        preemptible: 2
    }
    output {
        File m2_pon="${output_name}.pon.vcf"
    }
}