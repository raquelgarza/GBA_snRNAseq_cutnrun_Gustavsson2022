#!/bin/python
import truster
import os

path_parent = os.path.dirname(os.getcwd())
raw_path = "/projects/fs5/jakobssonlab/ASAP_Project_Raw_Data_pt2/"
lunarc = "config_files/lunarc_config.json"
modules = "config_files/software_modules.json"

asap = truster.Experiment("asap_5prime", lunarc, modules)

asap.register_sample(sample_id = "DA478", sample_name = "DA478-ASAP13-Ctl-NP16-161-PFCTX-5endseq-Seq135-4", raw_path = os.path.join(raw_path, "seq134B_135_137_139/fastq/ASAP/DA478-ASAP13-Ctl-NP16-161-PFCTX-5endseq-Seq135-4/"))
asap.register_sample(sample_id = "DA480", sample_name = "DA480-ASAP15-Ctl-NP16-21-PFCTX-5endseq-Seq135-6", raw_path = os.path.join(raw_path, "seq134B_135_137_139/fastq/ASAP/DA480-ASAP15-Ctl-NP16-21-PFCTX-5endseq-Seq135-6/"))
asap.register_sample(sample_id = "DA488", sample_name = "DA488-ASAP16-Ctl-PT231-PFCTX-5endseq-Seq135-8", raw_path = os.path.join(raw_path, "seq134B_135_137_139/fastq/ASAP/DA488-ASAP16-Ctl-PT231-PFCTX-5endseq-Seq135-8/"))
asap.register_sample(sample_id = "Seq127_1", sample_name = "ASAP13_Ctl_NP16-161_PFCTX_5prim", raw_path = "/projects/fs5/jakobssonlab/CTG_JGJSeq127_128/210907_A00681_0456_BH3KL5DMXY/H3KL5DMXY/outs/fastq_path/10X/Seq127_1/")

quantification_dir = os.path.join(path_parent, "1_counts")
cellranger_index = "/projects/fs1/common/genome/lunarc/10Xindexes/cellranger/6.0/refdata-gex-GRCh38-2020-A/" 
gene_gtf = "/projects/fs3/raquelgg/annotations/hg38/gencode/v38/gencode.v38.annotation.gtf"
te_gtf = "/projects/fs3/raquelgg/annotations/hg38/repeatmasker/hg38_rmsk_TEtranscripts.gtf"

asap.quantify(cr_index = cellranger_index, outdir=quantification_dir, jobs=5, nuclei={"DA478" : True, "DA480" : True, "DA488" : True, "Seq127_1" : True})
# In case you had already run cellranger and just want to set the quantification directory to an existing one - comment the line above, uncomment the line below
#for sample_id in list(asap.samples.keys()):
#    asap.set_quantification_outdir(sample_id = sample_id, cellranger_outdir = os.path.join(quantification_dir, sample_id))

merged_dir = os.path.join(path_parent, "3_combinedUMAP_perCluster") 
gene_gtf = "/projects/fs3/raquelgg/annotations/hg38/gencode/v38/gencode.v38.annotation.gtf"
te_gtf = "/projects/fs3/raquelgg/annotations/hg38/repeatmasker/hg38_rmsk_TEtranscripts.gtf"
star_index = "/projects/fs5/jakobssonlab/GRCh38.p13_gencode.v38_STAR/" 

clusters_dir = os.path.join(path_parent, "2_getClusters")
asap.get_clusters_all_samples(clusters_dir, perc_mitochondrial = 10, normalization_method = "CLR", max_size=2000, res = 0.5, jobs=5)
# In case you had already run get_clusters for all samples and just want to set the directory to an existing one - comment the line above, uncomment the line below
#asap.set_clusters_outdir(clusters_outdir = clusters_dir)

merged_dir = os.path.join(path_parent, "3_combinedUMAP_perCelltype") 
asap.merge_samples(merged_dir, "CLR", res = 0.5)
# In case you had already run merge_samples and just want to set the directory to an existing one - comment the line above, uncomment the line below
#asap.set_merge_samples_outdir(merged_dir)
merged_pipeline_dir = os.path.join(merged_dir, "clusterPipeline")

# Here I loaded the resulting object from merge_samples into Rstudio to annotate cell types and wrote a tsv file per cell type (per sample) with the cell barcodes.
# I also created the columns in metadata called "cellType" (with the annotated cell type) and "condition" which was set to "PFCTX_Ctl" for all cells
# I then came back to this script and ran the following to multimap per cell type:
asap.process_clusters(mode = "merged", outdir = merged_pipeline_dir, groups = {"PFCTX_Ctl" : ["DA478", "DA480", "DA488", "Seq127_1"]}, factor = ["cellType", "condition"], gene_gtf = gene_gtf, te_gtf = te_gtf, star_index = star_index, RAM = 48725506423, jobs=8, unique=False, tsv_to_bam = True, filter_UMIs = True, bam_to_fastq = True, concatenate_lanes = True, merge_clusters = True, map_cluster = True, TE_counts = False, normalize_TE_counts = False)
# In case you had already run process_clusters and just want to register the clusters - comment the line above, uncomment the line below
# asap.set_merge_clusters(outdir_merged_clusters = "../3_combinedUMAP_perCelltype/clusterPipeline/merged_cluster/", groups = {"PFCTX_Ctl" : ["DA478", "DA480", "DA488", "Seq127_1"]})




