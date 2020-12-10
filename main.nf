#!/usr/bin/env nextflow


import java.nio.file.Paths
import groovy.io.FileType

// nextflow version: 0.32.0

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                            FUNCTIONS                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

/* ~~~~~~~~~~~~~~~~~~~ */
def absolute_path(path) {
/* ~~~~~~~~~~~~~~~~~~~ */

	def f = new File(path)

	if ( ! f.isAbsolute() ) {
		return Paths.get(workflow.projectDir.toString(), path).toString()
	} else {
		return path
	}
}

/* ~~~~~~~~~~~~~~~~~~~~~ */
def check_file_path(path) {
/* ~~~~~~~~~~~~~~~~~~~~~ */

	// a java file object
	def f = new File(path)
	def abs_path = ""
	
	// get absolute path
	if ( ! f.isAbsolute() ) {
		abs_path = Paths.get(workflow.projectDir.toString(), path).toString()
	} else {
		abs_path = path
	}

	// is the file exists ?
	def file = new File(abs_path)
	assert file.exists() : "Error " + file + " does not exist."

	return file.getAbsolutePath()
}

/* ~~~~~~~~~~~~~~~~~~~~~ */
def is_single_end(design) {
/* ~~~~~~~~~~~~~~~~~~~~~ */

	// get header
	def header = ""
	new File(design).withReader{ header = it.readLine() }

	// check columns names
	def columns = header.split(",")
	def file_col = columns.findIndexOf{ it == "file" } >= 0 ? true : false 
	def file1_col = columns.findIndexOf{ it == "file1" } >= 0 ? true : false 
	def file2_col = columns.findIndexOf{ it == "file2" } >= 0 ? true : false

	// determine single or paired end
	if ( file1_col & file2_col & !file_col ) {
		return false
	} else if ( file_col & !file1_col & !file2_col ) {
		return true
	} else {
		println "Error the design file is not well formatted."
		System.exit(1)
	}
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
def get_rough_read_length(read_length) {
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	
	// three babs star indices have been built with these read length parameters
	def starIndexReadLengths = [50, 75, 100]
	
	// take the index with the closest read length to the experiment's
	def diffs = []
	starIndexReadLengths.each() { length ->
		diff = (length - read_length.toInteger()).abs()
		diffs.add(diff)
	}
	def index = diffs.findIndexValues() { i -> i == diffs.min() }[0]
	def rough_read_length = starIndexReadLengths[index.toInteger()]

	return rough_read_length
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        SESSION PARAMETERS                           -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

//publish_mode = "copy"
publish_mode = "symlink"
publish_overwrite = true

// modules
md_cellranger = "CellRanger/3.0.2-bcl2fastq-2.20.0"
md_r = "R/3.6.0-foss-2016b-BABS"

//////////////
// directories

input_dir = absolute_path("input")
output_dir = absolute_path(params.output_directory)

script_dir = Paths.get(input_dir, "scripts").toString()
r_script_dir = Paths.get(script_dir, "r").toString()
py_script_dir = Paths.get(script_dir, "py").toString()

// directories
//////////////


/////////
Channel//
/////////
	.fromPath( Paths.get(input_dir, "sequencing", "paths.txt").toString() )
	.splitCsv(header:false)
	.map{ it[0] }
	.filter{ it =~ /.*.fastq.gz$/ }
	.map{ new File(it).getBaseName().replaceAll("_.*", "") }
	.unique()
	.set{ samples }

/////////
Channel//
/////////
	.fromPath( Paths.get(input_dir, "sequencing", "paths.txt").toString() )
	.splitCsv(header:false)
	.map{ it[0] }
	.filter{ it =~ /.*.fastq.gz$/ }
	.map{ new File(it).getParentFile().getParent() }
	.unique()
	.map{ [ 1 , it ]}.groupTuple().map{ it[1] }.map{ it.join(",") }
	.set{ directories }

/////////
samples//
/////////
	.combine(directories)
	.set{ to_cellranger }

/* ~~~~~~~~~~~~~~~~~ */
def extractInfo(path) {
/* ~~~~~~~~~~~~~~~~~ */

	def f = new File(path).getBaseName()

	def limsid = f.replaceAll("_.*", "")
	def lane = f.replaceAll(".*_L(\\d+)_.*", "\$1")
	def sample = f.replaceAll(".*_S(\\d+)_.*", "\$1")

	def read_tag = ""
	if ( f =~ ".*_(R\\d)_.*" ) {
		read_tag = f.replaceAll(".*_(R\\d)_.*", "\$1")
	} else {
		read_tag = f.replaceAll(".*_(I\\d)_.*", "\$1")
	}

	def m = [
		"limsid": limsid,
		"read": read_tag,
		"sample": Integer.parseInt(sample),
		"lane": Integer.parseInt(lane)
	]

	return [ m , path ]
}

/* ~~~~~~~~~~~~~~~ */
def sortFastq(x, y) {
/* ~~~~~~~~~~~~~~~ */

	if ( x["sample"] == y["sample"] & x["lane"] == y["lane"] ) {
		return 0
	} else if ( x["sample"] == y["sample"] & x["lane"] > y["lane"] ) {
		return 1
	} else if ( x["sample"] == y["sample"] & x["lane"] < y["lane"] ) {
		return -1
	} else if ( x["sample"] > y["sample"] ) {
		return 1
	} else {
		return -1
	}
}

/////////
Channel//
/////////
	.fromPath( Paths.get(input_dir, "sequencing", "paths.txt").toString() )
	.splitCsv(header:false)
	.map{ it[0] }
	.filter{ it =~ /.*.fastq.gz$/ }
	.map{ extractInfo(it) }
	.map{ [ it[0]["limsid"] , it[0]["read"] , it ] }
	.groupTuple(by:[0,1])
	.map{[
		it[0],
		it[1],
		it[2].sort{ x, y -> sortFastq(x[0], y[0]) }
	]}
	.map{[ it[0], it[1], it[2].collect{ x -> x[1] } ]}
	.set{ fastq_to_merge }

/////////////////////////
def fileSorting(list) {//
/////////////////////////

	def R1 = ""
	def R2 = ""
	def I1 = ""

	list.each{
		if ( it =~ /.*_R1_.*/ ) {
			R1 = it
		} else if ( it =~ /.*_R2_.*/ ) {
			R2 = it
		} else {
			I1 = it
		}
	}

	return [ R1 , R2 , I1 ]
}

/////////
Channel//
/////////
	.fromPath( Paths.get(input_dir, "sequencing", "paths.txt").toString() )
	.splitCsv(header:false)
	.map{ it[0] }
	.filter{ it =~ /.*.fastq.gz$/ }
	.map{[
		new File(it).getBaseName().replaceAll("_[IR]\\d_\\d+.fastq", ""),
		new File(it).getName()
	]}
	.groupTuple()
	.map{ [ it[0] , fileSorting(it[1]).join("\t") ] }
	.set{ paired_files_rows }

///////////////////////////
process paired_file_row {//
///////////////////////////

	tag { name }

	executor "slurm"
	time "00:03:00"

	publishDir Paths.get(output_dir, "geo").toString(),
		mode: publish_mode,
		overwrite: publish_overwrite

	input:
		set val(name), val(row) from paired_files_rows
	
	output:
		file "${name}.txt" into paired_files
	
	shell:
		"""
		echo "${row}" > "${name}.txt"
		"""
}

////////////////////////////
process paired_files_tsv {//
////////////////////////////

	executor "slurm"
	time "00:03:00"

	publishDir Paths.get(output_dir, "geo").toString(),
		mode: publish_mode,
		overwrite: publish_overwrite

	input:
		file txt from paired_files.collect()
	
	output:
		file "paired_files.txt" into paired_files_tsv
	
	shell:
		"""
		cat $txt | sort > "paired_files.txt"
		"""
}

/////////
Channel//
/////////
	.fromPath( Paths.get(input_dir, "sequencing", "paths.txt").toString() )
	.splitCsv(header:false)
	.map{ it[0] }
	.filter{ it =~ /.*.fastq.gz$/ }
	.set{ fastqs }

//////////////////////
process copy_fastq {//
//////////////////////

	tag { filename }

	executor "slurm"
	time "00:03:00"

	publishDir Paths.get(output_dir, "geo").toString(),
		mode: publish_mode,
		overwrite: publish_overwrite

	input:
		val fastq from fastqs
	
	output:
		file "${filename}" into to_md5, to_read_length
	
	shell:

		filename = new File(fastq).getName()

		"""
		cp -vL "${fastq}" "${filename}"
		"""
}

/////////////////////
process md5_fastq {//
/////////////////////

	tag { filename }

	executor "slurm"
	time "00:03:00"

	input:
		file fastq from to_md5
	
	output:
		stdout into md5
	
	shell:

		filename = new File(fastq.toString()).getName()

		"""
		md5sum "${fastq}" | awk '{ print \$2 "\t" \$1 }'
		"""
}

///////////////////////
process read_length {//
///////////////////////

	tag { filename }

	executor "slurm"
	time "00:03:00"

	input:
		file fastq from to_read_length
	
	output:
		set val(filename), stdout into read_length
	
	shell:

		filename = new File(fastq.toString()).getName()

		"""
		zcat "${fastq}" | head | sed -n '4p' | wc -c
		"""
}


/////
md5//
/////
	.map{ it.replace("\n", "")
	.split("\t").collect() }
	.set{ md5 }


/////////////
read_length//
/////////////
	.map{ [ it[0] , it[1].replace("\n", "")+"bp" ] }
	.set{ read_length }

/////
md5//
/////
	.join(read_length)
	.map{[
		it[0],
		"fastq",
		it[1],
		"Illumina HiSeq 4000",
		it[2],
		"paired-end"
	]}
	.map{ [ it[0].replace(".fastq.gz", "") , it.join("\t") ] }
	.set{ raw_file_rows }

////////////////////////
process raw_file_row {//
////////////////////////

	tag { name }

	executor "slurm"
	time "00:03:00"

	publishDir Paths.get(output_dir, "geo").toString(),
		mode: publish_mode,
		overwrite: publish_overwrite

	input:
		set val(name), val(row) from raw_file_rows
	
	output:
		file "${name}.txt" into raw_files
	
	shell:
		"""
		echo "${row}" > "${name}.txt"
		"""
}

/////////////////////////
process raw_files_tsv {//
/////////////////////////

	executor "slurm"
	time "00:03:00"

	publishDir Paths.get(output_dir, "geo").toString(),
		mode: publish_mode,
		overwrite: publish_overwrite

	input:
		file txt from raw_files.collect()
	
	output:
		file "raw_files.txt" into raw_files_tsv
	
	shell:
		"""
		cat $txt | sort > "raw_files.txt"
		"""
}


//////////////////////
process cellranger {//
//////////////////////

	tag { sample }

	executor "slurm"
	time "20h"
	cpus 16
	module md_cellranger

	publishDir Paths.get(output_dir, "cellranger").toString(),
		mode: "copy",
		overwrite: publish_overwrite

	input:
		set val(sample), val(folders) from to_cellranger

	output:
		set val(sample), file("*") into cellranger
		set \
			val(sample),
			val("Cell Ranger barcodes file"),
			file("${sample}/outs/raw_feature_bc_matrix/barcodes.tsv.gz") \
			into barcodes
		set \
			val(sample),
			val("Cell Ranger features file"),
			file("${sample}/outs/raw_feature_bc_matrix/features.tsv.gz") \
			into features
		set \
			val(sample),
			val("Cell Ranger matrices file"),
			file("${sample}/outs/raw_feature_bc_matrix/matrix.mtx.gz") \
			into matrices

	shell:
		"""
		cellranger count \
			--id=${sample} \
			--transcriptome=${params.transcriptome} \
			--fastqs=${folders} \
			--sample=${sample} \
			--localcores=${task.cpus}
		"""
}


////////////
//barcodes//
////////////
//	.concat(features)
//	.concat(matrices)
//	.println()

///////////////////////
process copy_matrix {//
///////////////////////

	// it's because extracting that directory was rerunning cellranger, and i
	// didn't that happen

	tag { sample }

	executor "local"

	input:
		set val(sample), file(cellranger_dir) from cellranger

	output:
		set val(sample), file("*.matrix") into cellranger_matrix

	shell:
		"""
		cp -Rv ${cellranger_dir}/outs/raw_feature_bc_matrix ./${sample}.matrix
		"""
}

///////////////////
cellranger_matrix//
///////////////////
	.map{ [it[0], it[0] == "MCK332A1" ? "IgG control" : "anti-Skint1", it[1] ]}
	.set{ cellranger_matrix }

////////////////////////////////
process create_seurat_object {//
////////////////////////////////

	// it's because extracting that directory was rerunning cellranger, and i
	// didn't that happen

	tag { sample }
	
	//executor "slurm"
	executor "local"
	cpus 1
	time "30m"

	module md_r

	publishDir Paths.get(output_dir, "cellranger", "rda").toString(),
		mode: publish_mode,
		overwrite: publish_overwrite

	input:
		set \
			val(sample),
			val(condition),
			file(cellranger_dir) from cellranger_matrix

	output:
		set val(sample), file("*.rda") \
			into \
				seurat_object,
				seurat_object_matrix

	shell:
		"""
		Rscript ${r_script_dir}/seurat/create_object.r \
			--directory ${cellranger_dir} \
			--project ${params.project} \
			--sample ${sample} \
			--condition "${condition}" \
			--min-cells 3 \
			--min-features 200 \
			--rda-out ${sample}.rda
		"""
}

//////////////////////
process raw_matrix {//
//////////////////////

	tag { sample }
	
	//executor "slurm"
	executor "local"
	cpus 1
	time "30m"

	module md_r

	publishDir Paths.get(output_dir, "cellranger", "matrix").toString(),
		mode: publish_mode,
		overwrite: publish_overwrite

	input:
		set val(sample), file(rda) from seurat_object_matrix

	output:
		set val(sample), file("*.csv") into seurat_matrix

		//Rscript ${r_script_dir}/seurat/export_matrix.r ${rda} ${sample}.csv
	shell:
		"""
		Rscript ${r_script_dir}/seurat/export_matrix.r ${rda} ${sample}.csv
		"""
}

