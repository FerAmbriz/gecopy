# gecopy
gecopy is a pipeline flow that coordinates three well-known programs to analyze copy number variation and generate a consensus of reported variants. The program utilizes three different Docker containers, which are coordinated by a script for seamless user experience. By leveraging these containers, the program streamlines the process of automating the analysis and reporting of copy number variations. It compares the results obtained from two out of the three programs, ensuring accurate and reliable consensus reports. This automation framework significantly enhances user convenience, as the script facilitates the seamless integration and execution of the three containers, making the entire process more user-friendly.

## Docker images
The copy number variant callers were individually containerized using docker for greater portability.
* [CNVkit](https://hub.docker.com/r/ambrizbiotech/cnvkit)
* [ExomeDepth](https://hub.docker.com/r/ambrizbiotech/exomedepth)
* [Panelcn_mops](https://hub.docker.com/r/ambrizbiotech/panelcnmops)

## Installation
For the installation, the only thing necessary is to have docker installed. Once installed, the automount script of the different docker containers is moved to the PATH
```
git clone https://github.com/FerAmbriz/gecopy.git
cd gecopy/scr
sudo mv gecopy_docker /usr/bin/
```
## Usage
The use consists of a single command with different flags for the assignment of each input.
```
Usage
	gecopy_docker [options]
Options
	-i --input		Folder with BAM files
	-o --output		Output folder location
	-r --ref		Reference genome file
	-n --normal		Folder with BAM files of normals
	-c --bed-cnvkit		File with regions of interest
	-e --bed-exomedepth 	File with regions of interest 
	-p --bed-panelcnmops	File with regions of interest
	-m --mappable 		File with mappable regions
	-f --refFlat 		File refFlat
```
## Output
Finally, the output will be arranged individually for each subfolder created in the output directory with the name of each copy number variation program used.
