# 1000genome-sequential

### Requirements:
- Python 3.6 or later
- Mapplotlib: `pip install matplotlib`

### Instructions on how to run the sequential workflow:
```
$ git clone https://github.com/rafaelfsilva/1000genome-sequential
$ cd 1000genome-sequential
$ export PATH=$PATH:$PWD/bin
$ ./prepare_input.sh
$ ./1000genome-workflow-spec.py
```

### Instructions to run a single program of the workflow:

The following instructions assume at least a successful run of the workflows has been completed.

From the `1000genome-sequential` folder run the following:

```
$ cd data/20130502
$ individuals.py ALL.chr1.250000.vcf 1 1 1001 3000
```
