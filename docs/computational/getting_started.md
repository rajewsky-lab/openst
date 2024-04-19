# Getting started

After folowing the [experimental protocol](../experimental/getting_started.md), we provide the [`openst`](https://pypi.org/project/openst/)
python package for transforming the raw sequencing data into objects that can be used for spatial, single-cell
analysis, in four steps:

1. [Preprocessing of sequencing](preprocessing_sequencing.md)
2. [Pairwise alignment](pairwise_alignment.md): the spatial coordinates of transcriptomics data are aligned
    to tissue imaging.
3. [Segmentation and single-cell quantification](generate_expression_matrix.md): transcriptomic data
    are aggregated into single cells using the information from cell segmentation of tissue images.
4. [3D reconstruction](threed_reconstruction.md) of tissue imaging and transcriptome from serial sections.
   *We provide tutorials for interactive visualization of 3D data.*

## Installation

### with pip, <small>recommended</small>

The computational tools of the Open-ST workflow are published as a [Python package]
and can be installed with `pip`, ideally by using a [virtual environment].
Open up a terminal and install `openst` with:

``` sh
pip install openst
```

!!! tip
    If you don't have prior experience with Python, we recommend reading
    [Using Python's pip to Manage Your Projects' Dependencies], which is a really
    good introduction on the mechanics of Python package management and helps you
    troubleshoot if you run into errors.

  [Python package]: https://pypi.org/project/openst/
  [virtual environment]: https://realpython.com/what-is-pip/#using-pip-in-a-python-virtual-environment
  [Markdown]: https://python-markdown.github.io/
  [Pygments]: https://pygments.org/
  [Python Markdown Extensions]: https://facelessuser.github.io/pymdown-extensions/
  [Using Python's pip to Manage Your Projects' Dependencies]: https://realpython.com/what-is-pip/

!!! warning
    If you have a Mac with Apple Silicon (M1 or later), please install `openst` on a Rosetta environment 
    to ensure full compatibility (i.e., `openst manual_pairwise_aligner` might not work otherwise). 
    Also, the version of scikit-image pinned in the requirements might not be available for Apple Silicon. 
    Thus, create an environment and install dependencies as follows (assuming you have installed `conda`):

    ```bash
    CONDA_SUBDIR=osx-64 conda create -n openst python=3.11
    conda install scikit-image==0.19.3
    pip install openst
    ```
    
    When running `openst manual_pairwise_aligner` for the first time, startup time will be longer than
    usual. This is the expected behavior with osx-64 binaries (Rosetta). Also, you can install the optional
    `napari` (used in `openst preview`):

    ```bash
    pip install napari
    ```

### from git

`openst` can be directly installed from the source [GitHub repository]:
```
git clone https://github.com/rajewsky-lab/openst.git
pip install -e openst
```

  [GitHub repository]: https://github.com/rajewsky-lab/openst

### with docker

The official [Docker image] is a great way to get up and running in a few
minutes, as it comes with all dependencies pre-installed. Open up a terminal
and pull the image with:

``` sh
docker pull rajewsky/openst
```

The `openst` executable is provided as an entry point and `serve` is the
default command.

You can run a terminal for `openst` by running
```sh
docker run -it openst
```

Optionally, you might want to enable X11 redirection to enable GUI support (e.g., for the `openst manual_pairwise_aligner` tool).
Follow these steps:

1. Install X Server on your host machine if not already installed.
2. Allow connections to your X Server by running the following command on your host:
```sh
xhost +
```
3. Run the Docker container with the following additional options:
```sh
docker run -it --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix openst
```

Now, you can execute PyQt5-based applications, and the GUI will be displayed on your host machine.

!!! note
    - Ensure that the X Server on your host allows connections (xhost +) before running the container,
    - Make sure the necessary dependencies are installed on your host machine for PyQt5 applications.
    - Remember to close the X Server connections after using the container:
      ```sh
      xhost -
      ```

    [Docker image]: https://hub.docker.com/r/rajewsky/openst/