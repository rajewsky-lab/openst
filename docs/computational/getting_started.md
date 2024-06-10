# Getting started

After folowing the [experimental protocol](../experimental/getting_started.md), we provide the [`openst`](https://pypi.org/project/openst/)
python package for transforming the raw sequencing data into objects that can be used for spatial, single-cell
analysis, in five steps:

1. [Preprocessing of Open-ST spatial barcodes (capture area)](preprocessing_capture_area.md)
2. [Preprocessing of Open-ST transcriptomic library](preprocessing_openst_library.md)
3. [Pairwise alignment](pairwise_alignment.md): the spatial coordinates of transcriptomics data are aligned
    to tissue imaging.
4. [Segmentation and single-cell quantification](generate_expression_matrix.md): transcriptomic data
    are aggregated into single cells using the information from cell segmentation of tissue images.
5. [3D reconstruction](threed_reconstruction.md) of tissue imaging and transcriptome from serial sections.
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

!!! warning "Running on Apple Silicon-based Macs"
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

??? warning "Installed on Linux, accessing via SSH from Windows or macOS"
    If you install `openst` on Linux, and then use SSH from Windows or macOS to run it, you will need to have
    X11 redirection for the GUI-based components (i.e., `manual_pairwise_aligner` and `preview`).

    For example, you have to run SSH as:
    ```bash
    ssh -X user@server
    ```

    If you are using macOS, please download and install [X-quartz](https://www.xquartz.org).

    If you are using Windows, we recommend using [Tabby](https://tabby.sh) as the terminal, then 
    [VcXsrv](https://sourceforge.net/projects/vcxsrv/) as the X Server. You will need to setup a new SSH Profile with 
    X11 forwarding in Tabby.

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
and pull the image:

``` sh
docker pull rajewsky/openst
```

#### Get a `bash` terminal
In a terminal:
```sh
docker run -it --rm --entrypoint bash openst
```

#### Attach a local folder to `docker`:
In a terminal:
```bash
docker run -it --rm -v local_folder:/app/docker_folder --entrypoint bash openst
```

Make sure to replace `local_folder` (and optionally `docker_folder`).


#### Support GUI (`openst manual_pairwise_aligner` or `napari`)
In a terminal:
```bash
docker run -it --rm -p 9876:9876 openst
```

Then open [https://localhost:9876](https://localhost:9876) in a browser.
The video below shows an example launching `openst manual_pairwise_aligner` on
the browser (we recommend Google Chrome or similar).

<video loop autoplay muted playsinline width="800">
<source src='../../static/video/openst_docker_gui.webm' type="video/webm">
</video>

!!! warning "Running on Apple Silicon-based Macs"
    If you have a Mac with Apple Silicon (M1 or later), you need to configure `docker run` as:

    ```bash
    docker run --platform linux/amd64 # .. rest of the command
    ```
    This will use emulation of amd64 binaries, so performance might be slower than native.
    Please consider running on a Linux-based workstation, if possible, with a CUDA-supported
    GPU.

[Docker image]: https://hub.docker.com/r/rajewsky/openst/