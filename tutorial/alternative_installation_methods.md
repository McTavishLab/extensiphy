# Alternative Extensiphy Installation Methods

Extensiphy may be installed by other methods than just the Docker installation method described in the [main README](https://github.com/McTavishLab/extensiphy/blob/main/README.md#building-and-testing-your-own-extensiphy-docker-image). Below are two alternative methods for installing Extensiphy and its dependency programs.


## Installing dependencies with Anaconda

The other fast way to install the dependency programs of Extensiphy is to use the Conda package manager.
The Conda package manager is excellent because it handles installing dependency programs very well.
The steps for installing the Extensiphy dependencies are pretty straight forward so lets walk through them.

1. Clone the Extensiphy repository (if you havent done this already) We need the `install.yml` file.

```bash
$ git clone https://github.com/McTavishLab/extensiphy.git
$ cd extensiphy
```

2. Install Conda using the [miniconda](https://docs.conda.io/en/latest/miniconda.html) installer.

3. Once you've cloned the Extensiphy repository and have Conda installed, create an Extensiphy environment by running this command.

```bash
conda env create -f install.yml
```

4. Activate your environment.

```bash
$ conda activate extensiphy_env
```

5. Finally, you can test your installation by running the following command.

```bash
$ ./extensiphy.sh -a ./testdata/combo.fas -d ./testdata
```

Once Extensiphy has finished running on the test data, you should see a line saying:
```bash
Alignment file is: [path/to]/EP_output/outputs/extended.aln
```
Congratulations! You're install of Extensiphy is complete and you are ready to run analyses.


## Installing Dependencies With `apt-get`
Some of the dependency programs for Extensiphy have `apt-get` installation methods for applicable operating systems. The install for dendropy is included because its an easy install with Python. Run the following commands to install these dependencies easily. You'll notice that these commands do not install all of Extensiphy's dependency programs. This is unfortunate but you'll have to install the remaineder of the dependencies by hand.

```
apt-get install raxml
apt-get install seqtk
apt-get install samtools
apt-get install bcftools
pip install dendropy
```


## Installing Depedencies By Hand

Its completely possible to install all of the dependencies for Extensiphy by hand.
If you know how to add programs to your PATH (and know what a PATH is), I probably don't need to explain how to do this.
However, to make sure this overview of installation instructions is complete, here is a brief description.
You will need to:
1. Download the dependency programs listed in the [dependencies](https://github.com/McTavishLab/extensiphy#dependencies) section of the Extensiphy README.
2. Unzip and make sure all of the necessary programs are executable.
3. Add all of the dependency programs to your computers PATH.


Thats it! Using Extensiphy is limited to Linux at the moment. Using Ubuntu will ensure the smoothest performance. If you want to use another distro or another flavor of Debian, you'll have to make sure you install analogous one-liners and all that. You have been warned.
