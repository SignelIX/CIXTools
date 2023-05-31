import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
requirements = []
with open("requirements.txt", "r") as fh2:
    for line in fh2:
        requirements.append(line.strip())
setuptools.setup(
    name="CIxTools",
    version="1.0.0",
    author="Eric Sigel, Citadel Discovery",
    author_email="eric@citadeldiscovery.io",
    description="Cheminformatics Tools for combinatorial enumerations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SignelIX/CIXTools",
    packages=setuptools.find_packages(),
    install_requires = requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
)