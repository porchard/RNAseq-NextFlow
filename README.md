# NextFlow pipeline for paired-end RNA-seq data

## Dependencies
Singularity (v. 3) and NextFlow (>= v. 20.10.0). Containers with the software for each step are pulled from the Sylabs cloud library (https://cloud.sylabs.io/library).

## Configuration
Paths to various generic files (STAR indices and GTF files) must be included in the nextflow.config file -- check that file and change paths accordingly. STAR index must be compatible with STAR v. 2.7.9a.

You'll also need to set the params.results variable -- either in the nextflow.config file itself, or on the command line when you run the pipeline ('--results /path/to/results').

Lastly, you'll need to include information about each RNA-seq library, including the genome to which it should be mapped, and the paths to the fastq files for each readgroup. Organize this information in a JSON file, as in library-config.json.

## Running
Once you have all of the above information, you can run the pipeline as follows (in this case, indicating the path to the results on the command line):

```bash
nextflow run -params-file library-config.json --results /path/to/results /path/to/main.nf
```
