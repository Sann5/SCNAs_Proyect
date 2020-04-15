Inside these directories are the .sif files and the node attribute file for each subtype and the healthy phenotype. 

Pipeline description:
	- sif_processing_pipeline.R: Filters and normalizes the top 133,170 interactions for each subtype:

		Input:
		1) subtype.sif.sorted

		Outputs:
		1) subtype_norm_133k_interactions.sif


############################################################ WARNING ########################################################################

The original raw (the input files for the sif_processing_pipeline) data files are not contained within these directories beacuse of their size. They can be found in the following directory of INMEGENS's remote server:

username@castillo.cluster.inemen.gob.mx:/home/labs/csbig/subtipos_mama_2018/networks/

	- her2.sif.sorted
	- basal.sif.sorted
	- luma.sif.sorted
	- lumb.sif.sorted
	- healthy.sif.sorted

