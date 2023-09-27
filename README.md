# cMF_data
This directory accumulates the examples of cMF calculation based on Mayhall group's cMF [link](https://github.com/nmayhall-vt/ClusterMeanField.jl) implemenattaion.
Perform `CMF` (Cluster Mean-Field) calculations. This is simply a variational optimization of both orbital and cluster state parameters, minimizing the energy of a single TPS. This was originally proposed by Scuseria and coworkers [link](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.085101).

## Installation with Conda
1. Download

	```bash
	git clone https://github.com/nmayhall-vt/ClusterMeanField.jl.git
	cd ClusterMeanField.jl/
	```


2. Create conda environment to install Julia and will hold the PySCF executable. Install Julia with conda makes sure the correct python version will be found when using PyCall. where `-tauto` let's Julia pick the max number of threads. Use `-t N` to select `N` manually. Removing defaults to 1 thread. 

	```bash
	conda create -n my_env
	conda activate my_env
	conda install python==3.7
	conda config --add channels conda-forge
	conda install -c pyscf pyscf
	conda install h5py==2.10.0
	conda install julia
	julia --project=./ -tauto  
	```

3. Build PyCall from Julia REPL

  	```julia
	using Pkg; Pkg.build("PyCall")
	```


4. Run examples
	```
	include("./filepath")
	```

## Installation with Virtual Environment 
#### (might have issues with building PyCall and directing to correct python version)
1. Download

	```bash
	git clone git@github.com:nmayhall-vt/ClusterMeanField.jl.git
	cd ClusterMeanField.jl/
	```


2. Create python virtual environment which will hold the PYSCF executable.   where `-tauto` let's Julia pick the max number of threads. Use `-t N` to select `N` manually. Removing defaults to 1 thread. 

	```bash
	virtualenv -p python3 venv
	source venv/bin/activate
	pip install -r requirements.txt
	export PYTHON=$(which python)
	julia --project=./ -tauto 
  	```


3. Build PyCall from Julia REPL

  	```julia
	using Pkg; Pkg.build("PyCall")
	```
4. Run examples
	```
	include("./filepath")
	```
