# NextFlow pipeline for paired-end RNA-seq data

## Dependencies
If you have Singularity installed, you can use the config provided here ('Singularity') to build a container with all the dependencies.

Otherwise, you'll need to have the following installed:
1. STAR
2. fastqc
3. samtools
4. QoRTs 

I've used this pipeline with NextFlow v. 19.04.1

## Configuration
Paths to various generic files (STAR indices and chromosome size files) must be included in the nextflow.config file -- check that file and change paths accordingly.

You'll also need to set the params.results variable -- either in the nextflow.config file itself, or on the command line when you run the pipeline ('--results /path/to/results').

Lastly, you'll need to include information about each RNA-seq library, including the genome to which it should be mapped, and the paths to the fastq files for each readgroup. Organize this information in a JSON file, as in library-config.json.

## Running
Once you have all of the above information, you can run the pipeline as follows (in this case, indicating the path to the results on the command line):

```bash
nextflow run -with-singularity /path/to/Singularity.simg -params-file library-config.json --results /path/to/results /path/to/main.nf
```
