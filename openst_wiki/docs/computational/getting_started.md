# Getting started

After folowing the experimental protocol, we provide the [`openst`](https://pypi.org/project/openst/)
python package for transforming the raw sequencing data into objects that can be used for spatial, single-cell
analysis.

More specifically, our pipeline consists of the following steps:

1. [Preprocessing of sequencing](computational/preprocessing_sequencing.md)
2. [Preprocessing of imaging](computational/preprocessing_imaging.md)
3. [Align image to transcriptome](computational/pairwise_alignment.md): the spatial coordinates of transcripts are aligned
    to the imaging modality.
4. [Generating a cell-by-gene matrix](computational/generate_expression_matrix.md): transcripts
    are quantified per cell using the segmentation information.
5. [3D reconstruction](computational/threed_reconstruction.md) of tissue imaging and transcriptome from serial sections.
   *We provide tutorials for interactive visualization of 3D data.*

If you're familiar with Python, you can install `openst` with [pip], the Python package manager.
If not, we recommend using [docker].

[pip]: #with-pip
[docker]: #with-docker

## Installation

### with pip <small>recommended</small>

The computational tools of the open-ST workflow are published as a [Python package]
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

### with docker

The official [Docker image] is a great way to get up and running in a few
minutes, as it comes with all dependencies pre-installed. Open up a terminal
and pull the image with:

``` sh
docker pull rajewsky-lab/openst
```

The `openst` executable is provided as an entry point and `serve` is the
default command.

    [Docker image]: [https://docker](https://hub.docker.com/r/rajewsky-lab/openst/)

### with git

`openst` can be directly used from [GitHub] by cloning the
repository into a subfolder of your project root which might be useful if you
want to use the very latest version:

```
git clone https://github.com/rajewsky-lab/openst.git
```

Next, install the theme and its dependencies with:

```
pip install -e openst
```

  [GitHub]: https://github.com/rajewsky-lab/openst