# Modification_Effect_Socioeconomic_smoking_Asthma
This is a reproducible Bayesian Network analysis for the modification effect of Socioeconomic status on the effect of smoking on asthma related outcomes
If any interest in the reproducing the results please contact the author of this repo for more directions
In order to set the environment for reproducing the results in a Docker environment, one need to do the following:
- Git clone this repository
- Add a directory under the name data and upload the given data inside it.
- Add a directory inside Intermediate under the name outputs and add the give r objects. Specifically, the mice and bn objects.
- pull the following Docker image ranibasna/smok_soci_miss_bn_sing_renv
- open the command line and run the following command
- docker run -it --rm -v /path_to_your_local_project_directory:/socioeconomics ranibasna/smok_soci_miss_bn_sing_renv

## If you want to use the rstudio and skip the docker environment follow the below steps, However, without the Docker iamge the results may change a bit and it will not produce a full reproducible workflow.
- Git clone this repository
- Add a directory under the name data and upload the given data inside it.
- Add a directory inside Intermediate under the name outputs and add the give r objects. Specifically, the mice and bn objects.
- install the renv package
- run renv::init()
- run targets::tar_make()
