import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CIxTools",
    version="1.0.0",
    author="Eric Sigel, Citadel Discovery",
    author_email="eric@citadeldiscovery.io",
    description="Cheminformatics Tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/<>/<>",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)