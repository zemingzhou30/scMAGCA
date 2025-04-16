import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scMAGCA",
    version="1.0.10",
    author="zhouzeming",
    author_email="657608841@qq.com",
    description="scMAGCA, a multi-omics adversarial graph convolutional autoencoder framework to extract interpretable latent embeddings that characterize scMultiomics data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zemingzhou30/scMAGCA",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
