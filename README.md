# gimmecpg-r

Usage: gimmecpg.R [options]


Options:
        -i INPUT, --input=INPUT
                Path to bed files

        -n NAME, --name=NAME
                Name of specific files. Multiple filenames can be separated with a comma

        -e EXCLUDE, --exclude=EXCLUDE
                Path to a list of CpG sites to exclude

        -o OUTPUT, --output=OUTPUT
                Path to output directory

        -r REF, --ref=REF
                Path to reference methylation file

        -c MINCOV, --minCov=MINCOV
                Minimum coverage to consider methylation site as present. Default = 10

        -d MAXDISTANCE, --maxDistance=MAXDISTANCE
                Maximum distance between missing site and each neighbour for the site to be imputed. Default = 1000 bp

        -k, --collapse
                Choose whether to merge methylation sites on opposite strands together. Default = True

        -x, --machineLearning
                Choose whether to use machine learning for imputation. Default = False (no machine learning)

        -t RUNTIME, --runTime=RUNTIME
                Time (seconds) to train model. Default = 3600s (2h)

        -m MAXMODELS, --maxModels=MAXMODELS
                Maximum number of models to train within the time specified under --runTime. Excludes Stacked Ensemble models

        -s, --streaming
                Choose if streaming is required (for files that exceed memory). Default = False

        -h, --help
                Show this help message and exit