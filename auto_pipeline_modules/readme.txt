Idea.

Break automation of ELBA into three independent modules:

    1. automated compilation (grid sweep across k-mer bounds)
    2. automated "execution"
        a. automated slurm script generation for large jobs
        b. automated execution for running during interactive job
    3. automated post-processing of output.

In order to make the post-processing of ELBA's output simple, make sure all output
from a single program execution is contained in a single directory (hopefully a directory
with a decent name explaining the specific read datasets, parameters, etc.)

