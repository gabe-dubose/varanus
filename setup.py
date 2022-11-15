from importlib.metadata import entry_points
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="varanus",
    version="2022.07",
    author="Gabe DuBose",
    author_email="gabe.dubose.sci@gmail.com",
    description="A package for annotating genomic variants.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gabe-dubose/varanus",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = ['pyfaidx'],
    scripts = ['bin/varanus']
)