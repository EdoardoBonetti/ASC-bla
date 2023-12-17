"""Setup file for TomBino package."""
from skbuild import setup
from setuptools import find_packages


with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

# setup(
#    name="idgenerator",
#    version="0.0.10",
#    description="An id generator that generated various types and lengths ids",
#    package_dir={"": "app"},
#    packages=find_packages(where="app"),
#    long_description=long_description,
#    long_description_content_type="text/markdown",
#    url="https://github.com/ArjanCodes/wip/package",
#    author="ArjanCodes",
#    author_email="arjan@arjancodes.com",
#    license="MIT",
#    classifiers=[
#        "License :: OSI Approved :: MIT License",
#        "Programming Language :: Python :: 3.10",
#        "Operating System :: OS Independent",
#    ],
#    install_requires=["bson >= 0.5.10"],
#    extras_require={
#        "dev": ["pytest>=7.0", "twine>=4.0.2"],
#    },
#    python_requires=">=3.10",
# )


test_deps = ["pytest", "numpy"]
docs = ["sphinx", "myst-nb", "pandocfilters"]
benchmarks = ["numpy", "tqdm", "pandas", "plotly", "matplotlib"]
all_deps = test_deps + docs


extras = {
    "test": test_deps,
    "docs": docs,
    "benchmark": benchmarks,
    "all": all_deps,
}

setup(
    name="TomBino",
    version="0.0.1",
    description="A simple basic linear algebra library",
    author="J. Schoeberl, E. Bonetti",
    author_email="bonettiedo@gmail.com",
    license="MIT",
    #package_dir={"TomBino": "TomBino"},
    packages=["TomBino"],
    long_description_content_type="text/markdown",
    # url="https://github.com/",
    # classifiers=[
    #    "License :: OSI Approved :: MIT License",
    #    "Programming Language :: Python :: 3.10",
    #    "Operating System :: OS Independent",
    # ],
    cmake_args=[],
    # install_requires=["bson >= 0.5.10"],
    # extras_require={
    #    "dev": ["pytest>=7.0", "twine>=4.0.2"],
    # },
    # python_requires=">=3.10",
    tests_require=test_deps,
    extras_require=extras,
)
