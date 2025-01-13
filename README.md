Data and Julia codes for `Baseline dynamical model for three-versus-one ball possession in football' by H Ishii, Y Takai, et al.

# Setting up environment
You can reproduce the Julia environment with the required packages using Docker (`Dockerfile`).
The packages are installed referring to `Project.toml` and `Manifest.toml` while building a docker container.
- With VS Code `Dev Containers' extension and Docker Desktop [recommended]
  - Simply download this repository
  - Open the directory as container, using `Dev Containers: Open Folder in Container...` command.
    - You may need to reload the window after container is built.
    - Relevant VS Code extensions will be installed in the workspace.
      (Your global environment should NOT be affected.)
  - You can run the codes inside the container.
- With Docker Desktop
  - Use `Dockerfile` to reproduce the environment.
  - You may need to figure out how to run `.ipynb` files inside the container.
   
# Included files
- `src` : Directory containing the custom Julia module for running simulations and analysing resulting data
  - `BallPossessionModel.jl` : Main file of `BallPossessionModel` module
  - `BPMSimulation.jl` : Codes for performing simulations
  - `BPMAnalysis.jl` : Codes for analysing resulting data
- `scripts` : Directory within which codes should be run
  - `generate-dataset.jl` : Code to generate the simulated dataset which was run on HPC
  - `sim-dataset` : Directory for input and output of `generate-dataset.jl`
    - `sim-dataset-header.csv` : Input for `generate-dataset.jl`
    - `sim-dataset-numpass-data.csv` : Generated dataset
  - `plot-emp-dataset.ipynb` : Generate Figure 1
  - `emp-dataset` : Input of `plot-emp-dataset.ipynb`
    - `numpass-area-top.csv` : Data for the high-level team
    - `numpass-area-2nd.csv` : Data for the lower-level team
  - `emp-histograms.pdf` : Figure 1
  - `plot-sim-dataset.ipynb` : Generate Supplementary Figures S1 and S2 and Figure 3
  - `sim-dataset-hist-ord_data.pdf` : Supplementary Figure S1
  - `sim-dataset-hist-pro_data.pdf` : Supplementary Figure S2
  - `sim-dataset-area-numpass-scatter.pdf` : Figure 3
  - `regression.ipynb` : Generate Figure 4, Supplementary Figure S3, and Supplementary Tables S1
  - `sim-dataset-ols_qua-parameffect.pdf` : Figure 4
  - `sim-dataset-ols_qua-interaction.pdf` : Supplementary Figure S3
  - `sim-dataset-regtab.tex` : Original `.tex` file for Supplementary Table S1 (manually modified in manuscript)
  - `linear-reg-table.txt` : Full output on linear model
  - `square-reg-table.txt` : Full output on square model w/o interaction terms
  - `quadratic-reg-table.txt` : Full output on quadratic model
  - `generate-animation.ipynb` : Generate Supplementary Videos S1, S2, S3, and S4
  - `S1_intercept-w_guide.mp4` : Supplementary Video S1
  - `S2_maxpass-wo_guide.mp4` : Supplementary Video S2
  - `S3_ballout-wo_guide.mp4` : Supplementary Video S3
  - `S4_intercept-wo_guide.mp4` : Supplementary Video S4
- `.devcontainer` : For building a container with VS Code `Dev Containers' extension
  - `.DockerFile` : Dockerfile
  - `devcontainer.json` : Relevant to `Dev Containers' extension
- `hpc-init-sif.sh` : Shell script to initialize an apptainer container (`.sif`) on HPC
- `hpc-gen-data.sh` : Shell script to run `generate-dataset.jl` within the apptainer container on HPC
- `Manifest.toml` : Relevant to Julia's package management
- `Project.toml` : Relevant to Julia's package management
- `README.md` : This file
- `LICENSE` : License file

Should you have any question, please get in touch with me (ISHII Hidemasa).
