import io
import os
from setuptools import find_packages, setup

def read(*paths, **kwargs):
    content = ""
    with io.open(
        os.path.join(os.path.dirname(__file__), *paths),
        encoding=kwargs.get("encoding", "utf8"),
    ) as open_file:
        content = open_file.read().strip()
    return content

def read_requirements(path):
    return [
        line.strip()
        for line in read(path).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


setup(
    name="openst",
    version="0.1.0",
    description="Computational tools of the open-ST pipeline",
    author="Daniel León-Periñán",
    packages=find_packages(exclude=["tests", ".github", "models", "examples"]),
    entry_points = {
        'console_scripts': ['openst=openst.__main__:run_openst'],
    },
    package_data={'': ['openst/metadata/templates/*.html', '*/templates/*.html', 'openst/preprocessing/imagej_macros/*.ijm', '*/preprocessing/imagej_macros/*.ijm']},
    include_package_data=True,
)