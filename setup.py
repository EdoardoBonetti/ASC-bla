from skbuild import setup


setup(
    name="ASCsoft",
    version="0.0.2",
    author="J. Schoeberl, E. Bonetti",
    license="MIT",
    packages=["ASCsoft"],
    cmake_args=["-Wall"], 
)
